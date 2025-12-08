#!/usr/env python

"""
Pypiper is a python module with two components:
1) the PipelineManager class, and
2) other toolkits (currently just NGSTk) with functions for more specific pipeline use-cases.
The PipelineManager class can be used to create a procedural pipeline in python.
"""

import atexit
import datetime
import errno
import glob
import os
import platform
import re
import shlex  # for splitting commands like a shell does
import signal
import subprocess
import sys
import time
import warnings

from collections.abc import Iterable
from hashlib import md5

import __main__
import logmuse
import pandas as _pd
import psutil
from attmap import AttMapEcho
from pipestat import PipestatError, PipestatManager
from yacman import load_yaml

from ._version import __version__
from .const import PROFILE_COLNAMES, DEFAULT_SAMPLE_NAME
from .exceptions import PipelineHalt, SubprocessError
from .flags import *
from .utils import (
    CHECKPOINT_SPECIFICATIONS,
    check_shell,
    checkpoint_filepath,
    clear_flags,
    default_pipeline_config,
    flag_name,
    get_proc_name,
    is_multi_target,
    logger_via_cli,
    make_lock_name,
    parse_cmd,
    pipeline_filepath,
    default_pipestat_output_schema,
    result_formatter_markdown,
)
from yacman import load_yaml

__all__ = ["PipelineManager"]


LOCK_PREFIX = "lock."
LOGFILE_SUFFIX = "_log.md"


class Unbuffered(object):
    def __init__(self, stream):
        self.stream = stream

    def write(self, data):
        self.stream.write(data)
        self.stream.flush()

    def writelines(self, datas):
        self.stream.writelines(datas)
        self.stream.flush()

    def __getattr__(self, attr):
        return getattr(self.stream, attr)


class PipelineManager(object):
    """
    Base class for instantiating a PipelineManager object,
    the main class of Pypiper.

    :param str name: Choose a name for your pipeline;
        it's used to name the output files, flags, etc.
    :param str outfolder: Folder in which to store the results.
    :param argparse.Namespace args: Optional args object from ArgumentParser;
        Pypiper will simply record these arguments from your script
    :param bool multi: Enables running multiple pipelines in one script
        or for interactive use. It simply disables the tee of the output,
        so you won't get output logged to a file.
    :param bool dirty: Overrides the pipeline's clean_add()
        manual parameters, to *never* clean up intermediate files automatically.
        Useful for debugging; all cleanup files are added to manual cleanup script.
    :param bool recover: Specify recover mode, to overwrite lock files.
        If pypiper encounters a locked target, it will ignore the lock and
        recompute this step. Useful to restart a failed pipeline.
    :param bool new_start: start over and run every command even if output exists
    :param bool force_follow: Force run all follow functions
        even if  the preceding command is not run. By default,
        following functions  are only run if the preceding command is run.
    :param int cores: number of processors to use, default 1
    :param str mem: amount of memory to use. Default units are megabytes unless
        specified using the suffix [K|M|G|T]."
    :param str config_file: path to pipeline configuration file, optional
    :param str output_parent: path to folder in which output folder will live
    :param bool overwrite_checkpoints: Whether to override the stage-skipping
        logic provided by the checkpointing system. This is useful if the
        calls to this manager's run() method will be coming from a class
        that implements pypiper.Pipeline, as such a class will handle
        checkpointing logic automatically, and will set this to True to
        protect from a case in which a restart begins upstream of a stage
        for which a checkpoint file already exists, but that depends on the
        upstream stage and thus should be rerun if it's "parent" is rerun.
    :param str pipestat_record_identifier: record_identifier to report results via pipestat
    :param str pipestat_schema: output schema used by pipestat to report results
    :param str pipestat_results_file: path to file backend for reporting results
    :param str pipestat_config_file: path to pipestat configuration file
    :param str pipestat_pipeline_type: Sample or Project level pipeline
    :param pipestat_result_formatter: function used to style reported results, defaults to result_formatter_markdown
    :raise TypeError: if start or stop point(s) are provided both directly and
        via args namespace, or if both stopping types (exclusive/prospective
        and inclusive/retrospective) are provided.
    """

    # TODO: add pipestat-related args docstrings

    def __init__(
        self,
        name,
        outfolder,
        version=None,
        args=None,
        multi=False,
        dirty=False,
        recover=False,
        new_start=False,
        force_follow=False,
        cores=1,
        mem="1000M",
        config_file=None,
        output_parent=None,
        overwrite_checkpoints=False,
        logger_kwargs=None,
        pipestat_record_identifier=None,
        pipestat_schema=None,
        pipestat_results_file=None,
        pipestat_config=None,
        pipestat_pipeline_type=None,
        pipestat_result_formatter=None,
        **kwargs,
    ):
        # Params defines the set of options that could be updated via
        # command line args to a pipeline run, that can be forwarded
        # to Pypiper. If any pypiper arguments are passed
        # (via add_pypiper_args()), these will override the constructor
        # defaults for these arguments.

        # Establish default params
        params = {
            "dirty": dirty,
            "recover": recover,
            "new_start": new_start,
            "force_follow": force_follow,
            "config_file": config_file,
            "output_parent": output_parent,
            "cores": cores,
            "mem": mem,
            "testmode": False,
        }

        # Transform the command-line namespace into a Mapping.
        args_dict = vars(args) if args else dict()

        # Parse and store stage specifications that can determine pipeline
        # start and/or stop point.
        # First, add such specifications to the command-line namespace,
        # favoring the command-line spec if both are present.
        for cp_spec in set(CHECKPOINT_SPECIFICATIONS) & set(kwargs.keys()):
            args_dict.setdefault(cp_spec, kwargs[cp_spec])
        # Then, ensure that we set each such specification on this manager
        # so that we're guaranteed safe attribute access. If it's present,
        # remove the specification from the namespace that will be used to
        # update this manager's parameters Mapping.
        for optname in CHECKPOINT_SPECIFICATIONS:
            checkpoint = args_dict.pop(optname, None)
            setattr(self, optname, checkpoint)
        if self.stop_before and self.stop_after:
            raise TypeError(
                "Cannot specify both pre-stop and post-stop; "
                "got '{}' and '{}'".format(self.stop_before, self.stop_after)
            )

        # Update this manager's parameters with non-checkpoint-related
        # command-line parameterization.
        params.update(args_dict)

        # If no starting point was specified, assume that the pipeline's
        # execution is to begin right away and set the internal flag so that
        # run() is let loose to execute instructions given.
        self._active = not self.start_point

        # Pipeline-level variables to track global state and pipeline stats
        # Pipeline settings
        self.name = name
        self.tee = None
        self.overwrite_locks = params["recover"]
        self.new_start = params["new_start"]
        self.force_follow = params["force_follow"]
        self.dirty = params["dirty"]
        self.cores = params["cores"]
        self.output_parent = params["output_parent"]
        self.testmode = params["testmode"]

        # Establish the log file to check safety with logging keyword arguments.
        # Establish the output folder since it's required for the log file.
        self.outfolder = os.path.join(outfolder, "")  # trailing slash
        self.pipeline_log_file = pipeline_filepath(self, suffix=LOGFILE_SUFFIX)

        # Set up logger
        logger_kwargs = logger_kwargs or {}
        if logger_kwargs.get("logfile") == self.pipeline_log_file:
            raise ValueError(
                f"The logfile given for the pipeline manager's logger matches that which will be used by the manager itself: {self.pipeline_log_file}"
            )
        default_logname = ".".join([__name__, self.__class__.__name__, self.name])
        self._logger = None
        if args:
            logger_builder_method = "logger_via_cli"
            try:
                self._logger = logger_via_cli(args, **logger_kwargs)
            except logmuse.est.AbsentOptionException as e:
                # Defer logger construction to init_logger.
                self.debug(f"logger_via_cli failed: {e}")
        if self._logger is None:
            logger_builder_method = "init_logger"
            # covers cases of bool(args) being False, or failure of logger_via_cli.
            # strict is only for logger_via_cli.
            logger_kwargs = {k: v for k, v in logger_kwargs.items() if k != "strict"}
            try:
                name = logger_kwargs.pop("name")
            except KeyError:
                name = default_logname
            self._logger = logmuse.init_logger(name, **logger_kwargs)
        self.debug(f"Logger set with {logger_builder_method}")

        # Keep track of an ID for the number of processes attempted
        self.proc_count = 0

        # We use this memory to pass a memory limit to processes like java that
        # can take a memory limit, so they don't get killed by a SLURM (or other
        # cluster manager) overage. However, with java, the -Xmx argument can only
        # limit the *heap* space, not total memory use; so occasionally SLURM will
        # still kill these processes because total memory goes over the limit.
        # As a kind of hack, we'll set the java processes heap limit to 95% of the
        # total memory limit provided.
        # This will give a little breathing room for non-heap java memory use.

        if not params["mem"].endswith(("K", "M", "G", "T")):
            self.mem = params["mem"] + "M"
        else:
            # Assume the memory is in megabytes.
            self.mem = params["mem"]

        self.javamem = str(int(int(self.mem[:-1]) * 0.95)) + self.mem[-1:]

        self.container = None
        self.clean_initialized = False

        # Do some cores math for split processes
        # If a pipeline wants to run a process using half the cores, or 1/4 of the cores,
        # this can lead to complications if the number of cores is not evenly divisible.
        # Here we add a few variables so that pipelines can easily divide the cores evenly.
        # 50/50 split
        self.cores1of2a = int(self.cores) / 2 + int(self.cores) % 2
        self.cores1of2 = int(self.cores) / 2

        # 75/25 split
        self.cores1of4 = int(self.cores) / 4
        self.cores3of4 = int(self.cores) - int(self.cores1of4)

        self.cores1of8 = int(self.cores) / 8
        self.cores7of8 = int(self.cores) - int(self.cores1of8)

        self.pl_version = version
        # Set relative output_parent directory to absolute
        # not necessary after all. . .
        # if self.output_parent and not os.path.isabs(self.output_parent):
        #   self.output_parent = os.path.join(os.getcwd(), self.output_parent)

        # File paths:
        self.make_sure_path_exists(self.outfolder)
        self.pipeline_profile_file = pipeline_filepath(self, suffix="_profile.tsv")

        # Stats and figures are general and so lack the pipeline name.
        self.pipeline_stats_file = pipeline_filepath(self, filename="stats.yaml")

        # Record commands used and provide manual cleanup script.
        self.pipeline_commands_file = pipeline_filepath(self, suffix="_commands.sh")
        self.cleanup_file = pipeline_filepath(self, suffix="_cleanup.sh")

        # Pipeline status variables
        self.peak_memory = 0  # memory high water mark
        self.starttime = time.time()
        self.last_timestamp = self.starttime  # time of the last call to timestamp()

        self.locks = []
        self.running_procs = {}
        self.completed_procs = {}

        self.wait = True  # turn off for debugging

        # Initialize status and flags
        self.status = "initializing"
        # as part of the beginning of the pipeline, clear any flags set by
        # previous runs of this pipeline
        clear_flags(self)

        # In-memory holder for report_result
        self.stats_dict = {}

        # Result formatter to pass to pipestat
        self.pipestat_result_formatter = (
            pipestat_result_formatter or result_formatter_markdown
        )

        # Checkpoint-related parameters
        self.overwrite_checkpoints = overwrite_checkpoints or self.new_start
        self.halt_on_next = False
        self.prev_checkpoint = None
        self.curr_checkpoint = None

        # Pypiper can keep track of intermediate files to clean up at the end
        self.cleanup_list = []
        self.cleanup_list_conditional = []

        # Register handler functions to deal with interrupt and termination signals;
        # If received, we would then clean up properly (set pipeline status to FAIL, etc).
        signal.signal(signal.SIGINT, self._signal_int_handler)
        signal.signal(signal.SIGTERM, self._signal_term_handler)

        # pipestat setup
        self.pipestat_record_identifier = (
            pipestat_record_identifier or DEFAULT_SAMPLE_NAME
        )
        self.pipestat_pipeline_type = pipestat_pipeline_type or "sample"

        # don't force default pipestat_results_file value unless
        # pipestat config not provided
        if pipestat_config is None and pipestat_results_file is None:
            self.pipestat_results_file = self.pipeline_stats_file
        elif pipestat_results_file:
            self.pipestat_results_file = pipestat_results_file
            self.pipeline_stats_file = self.pipestat_results_file

        def _get_arg(args_dict, arg_name):
            """safely get argument from arg dict -- return None if doesn't exist"""
            return None if arg_name not in args_dict else args_dict[arg_name]

        self._pipestat_manager = PipestatManager(
            record_identifier=self.pipestat_record_identifier
            or _get_arg(args_dict, "pipestat_sample_name")
            or DEFAULT_SAMPLE_NAME,
            pipeline_name=self.name,
            schema_path=pipestat_schema
            or _get_arg(args_dict, "pipestat_schema")
            or default_pipestat_output_schema(sys.argv[0]),
            results_file_path=self.pipestat_results_file
            or _get_arg(args_dict, "pipestat_results_file")
            or self.pipeline_stats_file,
            config_file=pipestat_config or _get_arg(args_dict, "pipestat_config"),
            multi_pipelines=multi,
            pipeline_type=self.pipestat_pipeline_type,
        )

        self.start_pipeline(args, multi)

        # Handle config file if it exists

        # Read YAML config file
        # TODO: This section should become a function, so toolkits can use it
        # to locate a config file.
        config_to_load = None  # start with nothing

        if config_file:
            config_to_load = config_file
        else:
            cmdl_config_file = getattr(args, "config_file", None)
            if cmdl_config_file:
                if os.path.isabs(cmdl_config_file):
                    # Absolute custom config file specified
                    if os.path.isfile(cmdl_config_file):
                        config_to_load = cmdl_config_file
                    else:
                        self.debug("Can't find custom config file: " + cmdl_config_file)
                        pass
                else:
                    # Relative custom config file specified
                    # Set path to be relative to pipeline script
                    pipedir = os.path.dirname(sys.argv[0])
                    abs_config = os.path.join(pipedir, cmdl_config_file)
                    if os.path.isfile(abs_config):
                        config_to_load = abs_config
                    else:
                        self.debug("File: {}".format(__file__))
                        self.debug("Can't find custom config file: " + abs_config)
                        pass
                if config_to_load is not None:
                    pass
                    self.debug("\nUsing custom config file: {}".format(config_to_load))
            else:
                # No custom config file specified. Check for default
                default_config = default_pipeline_config(sys.argv[0])
                if os.path.isfile(default_config):
                    config_to_load = default_config
                    self.debug(
                        "Using default pipeline config file: {}".format(config_to_load)
                    )

        # Finally load the config we found.
        if config_to_load is not None:
            self.debug("\nLoading config file: {}\n".format(config_to_load))
            self.config = AttMapEcho(load_yaml(config_to_load))
        else:
            self.debug("No config file")
            self.config = None

    @property
    def pipestat(self):
        """
        `pipestat.PipestatManager` object to use for pipeline results reporting and status management

        Depending on the object configuration it can report to
        a YAML-formatted file or PostgreSQL database. Please refer to pipestat
        documentation for more details: http://pipestat.databio.org/

        :return pipestat.PipestatManager: object to use for results reporting
        """
        try:
            return getattr(self, "_pipestat_manager")
        except AttributeError:
            raise PipestatError(
                f"{PipestatManager.__name__} has not been configured for this pipeline run. "
                f"Provide an output schema to the {PipelineManager.__name__} object "
                f"in order to initialize it."
            )

    @property
    def _completed(self):
        """
        Is the managed pipeline in a completed state?

        :return bool: Whether the managed pipeline is in a completed state.
        """
        return (
            self.pipestat.get_status(self._pipestat_manager.record_identifier)
            == COMPLETE_FLAG
        )

    @property
    def _failed(self):
        """
        Is the managed pipeline in a failed state?

        :return bool: Whether the managed pipeline is in a failed state.
        """
        return (
            self.pipestat.get_status(self._pipestat_manager.record_identifier)
            == FAIL_FLAG
        )

    @property
    def halted(self):
        """
        Is the managed pipeline in a paused/halted state?
        :return bool: Whether the managed pipeline is in a paused/halted state.
        """
        return (
            self.pipestat.get_status(self._pipestat_manager.record_identifier)
            == PAUSE_FLAG
        )

    @property
    def _has_exit_status(self):
        """
        Has the managed pipeline been safely stopped?

        :return bool: Whether the managed pipeline's status indicates that it
            has been safely stopped.
        """
        return self._completed or self.halted or self._failed

    def _ignore_interrupts(self):
        """
        Ignore interrupt and termination signals. Used as a pre-execution
        function (preexec_fn) for subprocess.Popen calls that pypiper will
        control over (i.e., manually clean up).
        """
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGTERM, signal.SIG_IGN)

    def start_pipeline(self, args=None, multi=False):
        """
        Initialize setup. Do some setup, like tee output, print some diagnostics, create temp files.
        You provide only the output directory (used for pipeline stats, log, and status flag files).
        """
        # Perhaps this could all just be put into __init__, but I just kind of like the idea of a start function
        # self.make_sure_path_exists(self.outfolder)

        # By default, Pypiper will mirror every operation so it is displayed both
        # on sys.stdout **and** to a log file. Unfortunately, interactive python sessions
        # ruin this by interfering with stdout. So, for interactive mode, we do not enable
        # the tee subprocess, sending all output to screen only.
        # Starting multiple PipelineManagers in the same script has the same problem, and
        # must therefore be run in interactive_mode.

        interactive_mode = multi or not hasattr(__main__, "__file__")
        if interactive_mode:
            self.warning(
                "Warning: You're running an interactive python session. "
                "This works, but pypiper cannot tee the output, so results "
                "are only logged to screen."
            )
        else:
            sys.stdout = Unbuffered(sys.stdout)
            # sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)  # Unbuffer output

            # The tee subprocess must be instructed to ignore TERM and INT signals;
            # Instead, I will clean up this process in the signal handler functions.
            # This is required because otherwise, if pypiper receives a TERM or INT,
            # the tee will be automatically terminated by python before I have a chance to
            # print some final output (for example, about when the process stopped),
            # and so those things don't end up in the log files because the tee
            # subprocess is dead. Instead, I will handle the killing of the tee process
            # manually (in the exit handler).

            # a for append to file

            tee = subprocess.Popen(
                ["tee", "-a", self.pipeline_log_file],
                stdin=subprocess.PIPE,
                preexec_fn=self._ignore_interrupts,
            )

            # If the pipeline is terminated with SIGTERM/SIGINT,
            # make sure we kill this spawned tee subprocess as well.
            # atexit.register(self._kill_child_process, tee.pid, proc_name="tee")
            os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
            os.dup2(tee.stdin.fileno(), sys.stderr.fileno())

            self.tee = tee

        # For some reason, this exit handler function MUST be registered after
        # the one that kills the tee process.
        atexit.register(self._exit_handler)

        # A future possibility to avoid this tee, is to use a Tee class; this works for anything printed here
        # by pypiper, but can't tee the subprocess output. For this, it would require using threading to
        # simultaneously capture and display subprocess output. I shelve this for now and stick with the tee option.
        # sys.stdout = Tee(self.pipeline_log_file)

        # Record the git version of the pipeline and pypiper used. This gets (if it is in a git repo):
        # dir: the directory where the code is stored
        # hash: the commit id of the last commit in this repo
        # date: the date of the last commit in this repo
        # diff: a summary of any differences in the current (run) version vs. the committed version

        # Wrapped in try blocks so that the code will not fail if the pipeline or pypiper are not git repositories
        gitvars = {}
        try:
            # pypiper dir
            ppd = os.path.dirname(os.path.realpath(__file__))
            gitvars["pypiper_dir"] = ppd
            gitvars["pypiper_hash"] = (
                subprocess.check_output(
                    "cd " + ppd + "; git rev-parse --verify HEAD 2>/dev/null",
                    shell=True,
                )
                .decode()
                .strip()
            )
            gitvars["pypiper_date"] = (
                subprocess.check_output(
                    "cd " + ppd + "; git show -s --format=%ai HEAD 2>/dev/null",
                    shell=True,
                )
                .decode()
                .strip()
            )
            gitvars["pypiper_diff"] = (
                subprocess.check_output(
                    "cd " + ppd + "; git diff --shortstat HEAD 2>/dev/null", shell=True
                )
                .decode()
                .strip()
            )
            gitvars["pypiper_branch"] = (
                subprocess.check_output(
                    "cd " + ppd + "; git branch | grep '*' 2>/dev/null", shell=True
                )
                .decode()
                .strip()
            )
        except Exception:
            pass
        try:
            # pipeline dir
            pld = os.path.dirname(os.path.realpath(sys.argv[0]))
            gitvars["pipe_dir"] = pld
            gitvars["pipe_hash"] = (
                subprocess.check_output(
                    "cd " + pld + "; git rev-parse --verify HEAD 2>/dev/null",
                    shell=True,
                )
                .decode()
                .strip()
            )
            gitvars["pipe_date"] = (
                subprocess.check_output(
                    "cd " + pld + "; git show -s --format=%ai HEAD 2>/dev/null",
                    shell=True,
                )
                .decode()
                .strip()
            )
            gitvars["pipe_diff"] = (
                subprocess.check_output(
                    "cd " + pld + "; git diff --shortstat HEAD 2>/dev/null", shell=True
                )
                .decode()
                .strip()
            )
            gitvars["pipe_branch"] = (
                subprocess.check_output(
                    "cd " + pld + "; git branch | grep '*' 2>/dev/null", shell=True
                )
                .decode()
                .strip()
            )
        except Exception:
            pass

        # Print out a header section in the pipeline log:
        # Wrap things in backticks to prevent markdown from interpreting underscores as emphasis.
        # print("----------------------------------------")
        def logfmt(key, value=None, padding=16):
            padded_key = key.rjust(padding)
            formatted_val = f"`{value}`" if value else ""
            return f"* {padded_key}: {formatted_val}"

        self.info("### Pipeline run code and environment:\n")
        self.info(logfmt("Command", str(" ".join(sys.argv))))
        self.info(logfmt("Compute host", platform.node()))
        self.info(logfmt("Working dir", os.getcwd()))
        self.info(logfmt("Outfolder", self.outfolder))
        self.info(logfmt("Log file", self.pipeline_log_file))
        self.timestamp(logfmt("Start time"))

        self.info("\n### Version log:\n")
        self.info(logfmt("Python version", platform.python_version()))
        try:
            self.info(logfmt("Pypiper dir", gitvars["pypiper_dir"].strip()))
            self.info(logfmt("Pypiper version", __version__))
            self.info(logfmt("Pypiper hash", gitvars["pypiper_hash"]))
            self.info(logfmt("Pypiper branch", gitvars["pypiper_branch"]))
            self.info(logfmt("Pypiper date", gitvars["pypiper_date"]))
            if gitvars["pypiper_diff"]:
                self.info(logfmt("Pypiper diff", gitvars["pypiper_diff"]))
        except KeyError:
            # It is ok if keys aren't set, it means pypiper isn't in a  git repo.
            pass

        try:
            self.info(logfmt("Pipeline dir", gitvars["pipe_dir"].strip()))
            self.info(logfmt("Pipeline version", self.pl_version))
            self.info(logfmt("Pipeline hash", gitvars["pipe_hash"]).strip())
            self.info(logfmt("Pipeline branch", gitvars["pipe_branch"]).strip())
            self.info(logfmt("Pipeline date", gitvars["pipe_date"]).strip())
            if gitvars["pipe_diff"] != "":
                self.info(logfmt("Pipeline diff", gitvars["pipe_diff"]).strip())
        except KeyError:
            # It is ok if keys aren't set, it means the pipeline isn't a git repo.
            pass

        # self.info all arguments (if any)
        self.info("\n### Arguments passed to pipeline:\n")
        for arg, val in sorted((vars(args) if args else dict()).items()):
            argtext = "`{}`".format(arg)
            valtext = "`{}`".format(val)
            self.info("* {}:  {}".format(argtext.rjust(20), valtext))

        self.info("\n### Initialized Pipestat Object:\n")
        results = self._pipestat_manager.__str__().split("\n")
        for i in results:
            self.info("* " + i)
        self.info("* Sample name: " + self.pipestat_record_identifier + "\n")
        self.info("\n----------------------------------------\n")
        self.status = "running"
        self.pipestat.set_status(
            record_identifier=self._pipestat_manager.record_identifier,
            status_identifier="running",
        )

        # Record the start in PIPE_profile and PIPE_commands output files so we
        # can trace which run they belong to

        with open(self.pipeline_commands_file, "a") as myfile:
            myfile.write(
                "# Pipeline started at "
                + time.strftime("%m-%d %H:%M:%S", time.localtime(self.starttime))
                + "\n\n"
            )

        with open(self.pipeline_profile_file, "a") as myfile:
            myfile.write(
                "# Pipeline started at "
                + time.strftime("%m-%d %H:%M:%S", time.localtime(self.starttime))
                + "\n\n"
                + "# "
                + "\t".join(PROFILE_COLNAMES)
                + "\n"
            )

    def _set_status_flag(self, status):
        """
        Configure state and files on disk to match current processing status.

        :param str status: Name of new status designation for pipeline.
        """

        # Remove previous status flag file.
        flag_file_path = self._flag_file_path()
        try:
            os.remove(flag_file_path)
        except:
            # Print message only if the failure to remove the status flag
            # is unexpected; there's no flag for initialization, so we
            # can't remove the file.
            if self.status != "initializing":
                self.debug("Could not remove flag file: '{}'".format(flag_file_path))
            pass

        # Set new status.
        prev_status = self.status
        self.status = status
        self.pipestat.set_status(
            record_identifier=self._pipestat_manager.record_identifier,
            status_identifier=status,
        )
        self.debug("\nChanged status from {} to {}.".format(prev_status, self.status))

    def _flag_file_path(self, status=None):
        """
        Create path to flag file based on indicated or current status.

        Internal variables used are the pipeline name and the designated
        pipeline output folder path.

        :param str status: flag file type to create, default to current status
        :return str: path to flag file of indicated or current status.
        """

        flag_file_name = "{}_{}_{}".format(
            self._pipestat_manager.pipeline_name,
            self.pipestat_record_identifier,
            flag_name(status or self.status),
        )
        return pipeline_filepath(self, filename=flag_file_name)

    ###################################
    # Process calling functions
    ###################################
    def run(
        self,
        cmd,
        target=None,
        lock_name=None,
        shell=None,
        nofail=False,
        clean=False,
        follow=None,
        container=None,
        default_return_code=0,
    ):
        """
        The primary workhorse function of PipelineManager, this runs a command.

        This is the command  execution function, which enforces
        race-free file-locking, enables restartability, and multiple pipelines
        can produce/use the same files. The function will wait for the file
        lock if it exists, and not produce new output (by default) if the
        target output file already exists. If the output is to be created,
        it will first create a lock file to prevent other calls to run
        (for example, in parallel pipelines) from touching the file while it
        is being created. It also records the memory of the process and
        provides some logging output.

        :param str | list[str] cmd: Shell command(s) to be run.
        :param str | Sequence[str] target: Output file(s) to produce, optional.
            If all target files exist, the command will not be run. If no target
            is given, a lock_name must be provided.
        :param str lock_name: Name of lock file. Optional.
        :param bool shell: If command requires should be run in its own shell.
            Optional. Default: None --will try to determine whether the
            command requires a shell.
        :param bool nofail: Whether the pipeline proceed past a nonzero return from
            a process, default False; nofail can be used to implement
            non-essential parts of the pipeline; if a 'nofail' command fails,
            the pipeline is free to continue execution.
        :param bool clean: True means the target file will be automatically added
            to an auto cleanup list. Optional.
        :param callable follow: Function to call after executing (each) command.
        :param str container: Name for Docker container in which to run commands.
        :param Any default_return_code: Return code to use, might be used to discriminate
            between runs that did not execute any commands and runs that did.
        :return int: Return code of process. If a list of commands is passed,
            this is the maximum of all return codes for all commands.
        """

        def _max_ret_code(codes_list):
            """
            Return the maximum of a list of return codes.

            :param list[int] code: List of return codes to compare.
            :return int: Maximum of list.
            """
            # filter out codes that are None
            codes_list = [code for code in codes_list if code is not None]
            # get the max of the remaining codes
            if codes_list:
                return max(codes_list)
            # if no codes are left, return None
            return

        # validate default return code
        if default_return_code is not None and not isinstance(default_return_code, int):
            raise TypeError("default_return_code must be an int or None")

        # If the pipeline's not been started, skip ahead.
        if not self._active:
            cmds = [cmd] if isinstance(cmd, str) else cmd
            cmds_text = [c if isinstance(c, str) else " ".join(c) for c in cmds]
            self.info(
                "Pipeline is inactive; skipping {} command(s):\n{}".format(
                    len(cmds), "\n".join(cmds_text)
                )
            )
            return default_return_code

        # Short-circuit if the checkpoint file exists and the manager's not
        # been configured to overwrite such files.
        if self.curr_checkpoint is not None:
            check_fpath = checkpoint_filepath(self.curr_checkpoint, self)
            if os.path.isfile(check_fpath) and not self.overwrite_checkpoints:
                self.info(
                    "Checkpoint file exists for '{}' ('{}'), and the {} has "
                    "been configured to not overwrite checkpoints; "
                    "skipping command '{}'".format(
                        self.curr_checkpoint, check_fpath, self.__class__.__name__, cmd
                    )
                )
                return default_return_code

        # TODO: consider making the logic such that locking isn't implied, or
        # TODO (cont.): that we can make it otherwise such that it's not
        # TODO (cont.): strictly necessary to provide target or lock_name.
        # The default lock name is based on the target name.
        # Therefore, a targetless command that you want
        # to lock must specify a lock_name manually.
        if target is None and lock_name is None:
            self.fail_pipeline(
                Exception("You must provide either a target or a lock_name.")
            )

        # Downstream code requires target to be a list, so convert if only
        # a single item was given
        if not is_multi_target(target) and target is not None:
            target = [target]

        # Downstream code requires a list of locks; convert
        if isinstance(lock_name, str):
            lock_name = [lock_name]

        # Default lock_name (if not provided) is based on the target file name,
        # but placed in the parent pipeline outfolder
        self.debug(
            "Lock_name {}; target '{}', outfolder '{}'".format(
                lock_name, target, self.outfolder
            )
        )
        lock_name = lock_name or make_lock_name(target, self.outfolder)
        lock_files = [self._make_lock_path(ln) for ln in lock_name]

        process_return_code = default_return_code
        local_maxmem = 0

        # Decide how to do follow-up.
        if not follow:
            call_follow = lambda: None
        elif not hasattr(follow, "__call__"):
            # Warn about non-callable argument to follow-up function.
            self.warning(
                "Follow-up function is not callable and won't be used: {}".format(
                    type(follow)
                )
            )
            call_follow = lambda: None
        else:
            # Wrap the follow-up function so that the log shows what's going on.
            # additionally, the in_follow attribute is set to enable proper command count handling
            def call_follow():
                self.debug("Follow:")
                self.in_follow = True
                follow()
                self.in_follow = False

        # The while=True loop here is unlikely to be triggered, and is just a
        # wrapper to prevent race conditions; the lock_file must be created by
        # the current loop. If not, we loop again and then re-do the tests.
        # The recover and newstart options inform the pipeline to run a command
        # in a scenario where it normally would not. We use these "local" flags
        # to allow us to report on the state of the pipeline in the first round
        # as normal, but then proceed on the next iteration through the outer
        # loop. The proceed_through_locks is a flag that is set if any lockfile
        # is found that needs to be recovered or overwritten. It instructs us to
        # ignore lock files on the next iteration.
        local_recover = False
        local_newstart = False
        proceed_through_locks = False

        while True:
            ##### Tests block
            # Base case: All targets exists and not set to overwrite targets break loop, don't run process.
            # os.path.exists returns True for either a file or directory; .isfile is file-only
            if (
                target is not None
                and all([os.path.exists(t) for t in target])
                and not any([os.path.isfile(l) for l in lock_files])
                and not local_newstart
            ):
                for tgt in target:
                    if os.path.exists(tgt):
                        self.info("Target exists: `" + tgt + "`  ")
                if self.new_start:
                    self.info("New start mode; run anyway.  ")
                    # Set the local_newstart flag so the command will run anyway.
                    # Doing this in here instead of outside the loop allows us
                    # to still report the target existence.
                    local_newstart = True
                    continue
                # Normally we don't run the follow, but if you want to force. . .
                if self.force_follow:
                    call_follow()
                # Increment process count
                increment_info_pattern = "Skipped command: `{}`\nCommand ID incremented by: `{}`. Current ID: `{}`\n"
                if isinstance(cmd, list):
                    for c in cmd:
                        count = len(parse_cmd(c, shell))
                        self.proc_count += count
                        self.debug(
                            increment_info_pattern.format(
                                str(c), count, self.proc_count
                            )
                        )
                else:
                    count = len(parse_cmd(cmd, shell))
                    self.proc_count += count
                    self.debug(
                        increment_info_pattern.format(str(cmd), count, self.proc_count)
                    )
                break  # Do not run command

            # Scenario 1: Lock file exists, but we're supposed to overwrite target; Run process.
            if not proceed_through_locks:
                for lock_file in lock_files:
                    recover_file = self._recoverfile_from_lockfile(lock_file)
                    if os.path.isfile(lock_file):
                        self.info("Found lock file: {}".format(lock_file))
                        if self.overwrite_locks:
                            self.info("Overwriting target...")
                            proceed_through_locks = True
                        elif os.path.isfile(recover_file):
                            self.info(
                                "Found dynamic recovery file ({}); "
                                "overwriting target...".format(recover_file)
                            )
                            # remove the lock file which will then be promptly re-created for the current run.
                            local_recover = True
                            proceed_through_locks = True
                            # the recovery flag is now spent; remove so we don't accidentally re-recover a failed job
                            os.remove(recover_file)
                        else:  # don't overwrite locks
                            self._wait_for_lock(lock_file)
                            # when it's done loop through again to try one more
                            # time (to see if the target exists now)
                            continue

            # If you get to this point, the target doesn't exist, and the lock_file doesn't exist
            # (or we should overwrite). create the lock (if you can)
            # Initialize lock in master lock list
            for lock_file in lock_files:
                self.locks.append(lock_file)
                if self.overwrite_locks or local_recover:
                    self._create_file(lock_file)
                else:
                    try:
                        self._create_file_racefree(lock_file)  # Create lock
                    except OSError as e:
                        if e.errno == errno.EEXIST:  # File already exists
                            self.info(
                                "Lock file created after test! Looping again: {}".format(
                                    lock_file
                                )
                            )

                            # Since a lock file was created by a different source,
                            # we need to reset this flag to re-check the locks.
                            proceed_through_locks = False
                            continue  # Go back to start

            ##### End tests block
            # If you make it past these tests, we should proceed to run the process.

            if target is not None:
                self.info(
                    "Target to produce: {}  ".format(
                        ",".join(["`" + x + "`" for x in target])
                    )
                )
            else:
                self.info("Targetless command, running...  ")

            if isinstance(cmd, list):  # Handle command lists
                for cmd_i in cmd:
                    list_ret, maxmem = self.callprint(
                        cmd_i, shell, lock_file, nofail, container
                    )
                    maxmem = max(maxmem) if isinstance(maxmem, Iterable) else maxmem
                    local_maxmem = max(local_maxmem, maxmem)
                    list_ret = (
                        _max_ret_code(list_ret)
                        if isinstance(list_ret, Iterable)
                        else list_ret
                    )
                    process_return_code = _max_ret_code([process_return_code, list_ret])

            else:  # Single command (most common)
                process_return_code, local_maxmem = self.callprint(
                    cmd, shell, lock_file, nofail, container
                )  # Run command
                if isinstance(process_return_code, list):
                    process_return_code = _max_ret_code(process_return_code)

            # For temporary files, you can specify a clean option to automatically
            # add them to the clean list, saving you a manual call to clean_add
            if target is not None and clean:
                for tgt in target:
                    self.clean_add(tgt)

            call_follow()
            for lock_file in lock_files:
                os.remove(lock_file)  # Remove lock file
                self.locks.remove(lock_file)

            # If you make it to the end of the while loop, you're done
            break

        return process_return_code

    def checkprint(self, cmd, shell=None, nofail=False):
        """
        Just like callprint, but checks output -- so you can get a variable
        in python corresponding to the return value of the command you call.
        This is equivalent to running subprocess.check_output()
        instead of subprocess.call().
        :param str | Iterable[str] cmd: Bash command(s) to be run.
        :param bool | str shell: If command requires should be run in its own shell. Optional.
            Default: "guess" -- `run()` will try to guess if the command should be run in a
            shell (based on the presence of a pipe (|) or redirect (>), To force a process to
            run as a direct subprocess, set `shell` to False; to force a shell, set True.
        :param bool nofail: Should the pipeline bail on a nonzero return from a process? Default: False
            Nofail can be used to implement non-essential parts of the pipeline; if these processes fail,
            they will not cause the pipeline to bail out.
        :return str: text output by the executed subprocess (check_output)
        """

        self._report_command(cmd)
        if self.testmode:
            return ""

        likely_shell = check_shell(cmd, shell)

        if shell is None:
            shell = likely_shell

        if not shell:
            if likely_shell:
                self.debug(
                    "Should this command run in a shell instead of directly in a subprocess?"
                )
            cmd = shlex.split(cmd)

        try:
            return subprocess.check_output(cmd, shell=shell).decode().strip()
        except Exception as e:
            self._triage_error(e, nofail)

    def _attend_process(self, proc, sleeptime):
        """
        Waits on a process for a given time to see if it finishes, returns True
        if it's still running after the given time or False as soon as it
        returns.

        :param psutil.Popen proc: Process object opened by psutil.Popen()
        :param float sleeptime: Time to wait
        :return bool: True if process is still running; otherwise false
        """
        # print("attend:{}".format(proc.pid))
        try:
            proc.wait(timeout=sleeptime)
        except psutil.TimeoutExpired:
            return True
        return False

    def callprint(self, cmd, shell=None, lock_file=None, nofail=False, container=None):
        """
        Prints the command, and then executes it, then prints the memory use and
        return code of the command.

        Uses python's subprocess.Popen() to execute the given command. The shell argument is simply
        passed along to Popen(). You should use shell=False (default) where possible, because this enables memory
        profiling. You should use shell=True if you require shell functions like redirects (>) or stars (*), but this
        will prevent the script from monitoring memory use. The pipes (|) will be used to split the command into
        subprocesses run within python, so the memory profiling is possible.

        cmd can also be a series (a dict object) of multiple commands, which will be run in succession.

        :param str | Iterable[str] cmd: Bash command(s) to be run.
        :param str lock_file: a lock file name
        :param bool nofail: Should the pipeline bail on a nonzero return from a process? Default: False
            Nofail can be used to implement non-essential parts of the pipeline; if these processes fail,
            they will not cause the pipeline to bail out.
        :param bool shell: if the command should be run it its own shell, default: None (will try
            to determine based on the command)
        :param container: Named Docker container in which to execute.
        :param container: str
        """
        # The Popen shell argument works like this:
        # if shell=False, then we format the command (with split()) to be a list of command and its arguments.
        # Split the command to use shell=False;
        # leave it together to use shell=True;

        def get_mem_child_sum(proc):
            try:
                # get children processes
                children = proc.children(recursive=True)
                # get RSS memory of each child proc and sum all
                mem_sum = proc.memory_info().rss
                if children:
                    mem_sum += sum([x.memory_info().rss for x in children])
                # return in gigs
                return mem_sum / 1e9
            except (psutil.NoSuchProcess, psutil.ZombieProcess) as e:
                self.warning(e)
                self.warning(
                    "Warning: couldn't add memory use for process: {}".format(proc.pid)
                )
                return 0

        def display_memory(memval):
            return None if memval < 0 else "{}GB".format(round(memval, 3))

        def make_hash(o):
            """
            Convert the object to string and hash it, return None in case of failure
            :param o: object of any type, in our case it is a dict
            :return str: hahsed string representation of the dict
            """
            try:
                hsh = md5(str(o).encode("utf-8")).hexdigest()[:10]
            except Exception as e:
                self.debug(
                    "Could not create hash for '{}', caught exception: {}".format(
                        str(o), e.__class__.__name__
                    )
                )
                hsh = None
            return hsh

        if container:
            cmd = "docker exec " + container + " " + cmd

        if self.testmode:
            self._report_command(cmd)
            return 0, 0

        self.debug("Command: {}".format(cmd))
        param_list = parse_cmd(cmd, shell)
        # cast all commands to str and concatenate for hashing
        conc_cmd = "".join([str(x["args"]) for x in param_list])
        self.debug("Hashed command '{}': {}".format(conc_cmd, make_hash(conc_cmd)))
        processes = []
        running_processes = []
        completed_processes = []
        start_time = time.time()
        for i in range(len(param_list)):
            running_processes.append(i)
            if i == 0:
                processes.append(psutil.Popen(preexec_fn=os.setsid, **param_list[i]))
            else:
                param_list[i]["stdin"] = processes[i - 1].stdout
                processes.append(psutil.Popen(preexec_fn=os.setsid, **param_list[i]))
            self.running_procs[processes[-1].pid] = {
                "proc_name": get_proc_name(param_list[i]["args"]),
                "start_time": start_time,
                "container": container,
                "p": processes[-1],
                "args_hash": make_hash(conc_cmd),
                "local_proc_id": self.process_counter(),
            }

        self._report_command(cmd, [x.pid for x in processes])
        # Capture the subprocess output in <pre> tags to make it format nicely
        # if the markdown log file is displayed as HTML.
        self.info("<pre>")

        local_maxmems = [-1] * len(running_processes)
        returncodes = [None] * len(running_processes)
        proc_wrapup_text = [None] * len(running_processes)

        if not self.wait:
            self.info("</pre>")
            ids = [x.pid for x in processes]
            self.debug("Not waiting for subprocesses: " + str(ids))
            return 0, -1

        def proc_wrapup(i):
            """
            :param i: internal ID number of the subprocess
            """
            returncode = processes[i].returncode
            current_pid = processes[i].pid

            info = "PID: {pid};\tCommand: {cmd};\tReturn code: {ret};\tMemory used: {mem}".format(
                pid=current_pid,
                cmd=self.running_procs[current_pid]["proc_name"],
                ret=processes[i].returncode,
                mem=display_memory(local_maxmems[i]),
            )

            # report process profile
            self._report_profile(
                self.running_procs[current_pid]["proc_name"],
                lock_file,
                time.time() - self.running_procs[current_pid]["start_time"],
                local_maxmems[i],
                current_pid,
                self.running_procs[current_pid]["args_hash"],
                self.running_procs[current_pid]["local_proc_id"],
            )

            # Remove this as a running subprocess
            self.running_procs[current_pid]["info"] = info
            self.running_procs[current_pid]["returncode"] = returncode
            self.completed_procs[current_pid] = self.running_procs[current_pid]
            del self.running_procs[current_pid]
            running_processes.remove(i)
            completed_processes.append(i)
            proc_wrapup_text[i] = info
            returncodes[i] = returncode
            return info

        sleeptime = 0.0001

        while running_processes:
            self.debug("running")
            for i in running_processes:
                local_maxmems[i] = max(
                    local_maxmems[i], (get_mem_child_sum(processes[i]))
                )
                self.peak_memory = max(self.peak_memory, local_maxmems[i])
                self.debug(processes[i])
                if not self._attend_process(processes[i], sleeptime):
                    proc_wrapup_text[i] = proc_wrapup(i)

            # the sleeptime is extremely short at the beginning and gets longer exponentially
            # (+ constant to prevent copious checks at the very beginning)
            # = more precise mem tracing for short processes
            sleeptime = min((sleeptime + 0.25) * 3, 60 / len(processes))

        # All jobs are done, print a final closing and job info
        info = (
            "Elapsed time: "
            + str(datetime.timedelta(seconds=self.time_elapsed(start_time)))
            + "."
        )
        info += " Running peak memory: {pipe}.".format(
            pipe=display_memory(self.peak_memory)
        )
        # if len(proc_wrapup_text) == 1:
        # info += " {}".format(proc_wrapup_text[0])

        for i in completed_processes:
            info += "  \n  {}".format(self.completed_procs[processes[i].pid]["info"])

        info += "\n"  # finish out the
        self.info("</pre>")
        self.info("Command completed. {info}".format(info=info))

        for rc in returncodes:
            if rc != 0:
                msg = (
                    "Subprocess returned nonzero result. Check above output for details"
                )
                self._triage_error(SubprocessError(msg), nofail)

        return [returncodes, local_maxmems]

    def process_counter(self):
        """
        Increments process counter with regard to the follow state:
        if currently executed command is a follow function of another one, the counter is not incremented.

        :return str | int: current counter state, a number if the counter has beed incremented or a number
        of the previous process plus "f" otherwise
        """
        try:
            if self.in_follow:
                return str(self.proc_count) + "f"
            else:
                self.proc_count += 1
                return self.proc_count
        except AttributeError:
            self.proc_count += 1
            return self.proc_count

    ###################################
    # Waiting functions
    ###################################

    def _wait_for_process(self, p, shell=False):
        """
        Debug function used in unit tests.

        :param p: A subprocess.Popen process.
        :param bool shell: If command requires should be run in its own shell. Optional. Default: False.
        """
        local_maxmem = -1
        sleeptime = 0.5
        while p.poll() is None:
            if not shell:
                local_maxmem = max(local_maxmem, self._memory_usage(p.pid) / 1e6)
                # print("int.maxmem (pid:" + str(p.pid) + ") " + str(local_maxmem))
            time.sleep(sleeptime)
            sleeptime = min(sleeptime + 5, 60)

        self.peak_memory = max(self.peak_memory, local_maxmem)

        del self.running_procs[p.pid]

        info = "Process " + str(p.pid) + " returned: (" + str(p.returncode) + ")."
        if not shell:
            info += " Peak memory: (Process: " + str(round(local_maxmem, 3)) + "GB;"
            info += " Pipeline: " + str(round(self.peak_memory, 3)) + "GB)\n"

        self.info(info + "\n")
        if p.returncode != 0:
            raise Exception("Process returned nonzero result.")
        return [p.returncode, local_maxmem]

    def _wait_for_lock(self, lock_file):
        """
        Just sleep until the lock_file does not exist or a lock_file-related dynamic recovery flag is spotted

        :param str lock_file: Lock file to wait upon.
        """
        sleeptime = 0.5
        first_message_flag = False
        long_message_flag = False
        dot_count = 0
        totaltime = 0
        recover_file = self._recoverfile_from_lockfile(lock_file)
        while os.path.isfile(lock_file):
            if first_message_flag is False:
                self.timestamp("Waiting for file lock: " + lock_file)
                self.warning(
                    "This indicates that another process may be executing this "
                    "command, or the pipeline was not properly shut down.  If the "
                    "pipeline was not properly shut down last time, "
                    "you should restart it in 'recover' mode (-R) to indicate that "
                    "this step should be restarted."
                )
                # self._set_status_flag(WAIT_FLAG)
                self.pipestat.set_status(
                    record_identifier=self._pipestat_manager.record_identifier,
                    status_identifier="waiting",
                )
                first_message_flag = True
            else:
                sys.stdout.write(".")
                dot_count = dot_count + 1
                if dot_count % 60 == 0:
                    self.info("")  # linefeed
            # prevents the issue of pypiper waiting for the lock file to be gone infinitely
            # in case the recovery flag is sticked by other pipeline when it's interrupted
            if os.path.isfile(recover_file):
                sys.stdout.write(" Dynamic recovery flag found")
                break
            time.sleep(sleeptime)
            totaltime += sleeptime
            sleeptime = min(sleeptime + 2.5, 60)
            if sleeptime > 3600 and not long_message_flag:
                long_message_flag = True

        if first_message_flag:
            self.timestamp("File unlocked.")
            # self._set_status_flag(RUN_FLAG)
            self.pipestat.set_status(
                record_identifier=self._pipestat_manager.record_identifier,
                status_identifier="running",
            )

    ###################################
    # Logging functions
    ###################################

    def debug(self, msg, *args, **kwargs):
        self._logger.debug(msg, *args, **kwargs)

    def info(self, msg, *args, **kwargs):
        self._logger.info(msg, *args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        self._logger.warning(msg, *args, **kwargs)

    def error(self, msg, *args, **kwargs):
        self._logger.error(msg, *args, **kwargs)

    def critical(self, msg, *args, **kwargs):
        self._logger.critical(msg, *args, **kwargs)

    def fatal(self, msg, *args, **kwargs):
        self._logger.fatal(msg, *args, **kwargs)

    def timestamp(self, message="", checkpoint=None, finished=False, raise_error=True):
        """
        Print message, time, and time elapsed, perhaps creating checkpoint.

        This prints your given message, along with the current time, and time
        elapsed since the previous timestamp() call.  If you specify a
        HEADING by beginning the message with "###", it surrounds the message
        with newlines for easier readability in the log file. If a checkpoint
        is designated, an empty file is created corresponding to the name
        given. Depending on how this manager's been configured, the value of
        the checkpoint, and whether this timestamp indicates initiation or
        completion of a group of pipeline steps, this call may stop the
        pipeline's execution.

        :param str message: Message to timestamp.
        :param str checkpoint: Name of checkpoint; this tends to be something
            that reflects the processing logic about to be or having just been
            completed. Provision of an argument to this parameter means that
            a checkpoint file will be created, facilitating arbitrary starting
            and stopping point for the pipeline as desired.
        :param bool finished: Whether this call represents the completion of a
            conceptual unit of a pipeline's processing
        :param raise_error: Whether to raise exception if
            checkpoint or current state indicates that a halt should occur.
        """

        # Halt if the manager's state has been set such that this call
        # should halt the pipeline.
        if self.halt_on_next:
            self.halt(checkpoint, finished, raise_error=raise_error)

        # Determine action to take with respect to halting if needed.
        if checkpoint:
            if finished:
                # Write the file.
                self._checkpoint(checkpoint)
                self.prev_checkpoint = checkpoint
                self.curr_checkpoint = None
            else:
                self.prev_checkpoint = self.curr_checkpoint
                self.curr_checkpoint = checkpoint
                self._checkpoint(self.prev_checkpoint)
            # Handle the two halting conditions.
            if (finished and checkpoint == self.stop_after) or (
                not finished and checkpoint == self.stop_before
            ):
                self.halt(checkpoint, finished, raise_error=raise_error)
            # Determine if we've started executing.
            elif checkpoint == self.start_point:
                self._active = True
            # If this is a prospective checkpoint, set the current checkpoint
            # accordingly and whether we should halt the pipeline on the
            # next timestamp call.
            if not finished and checkpoint == self.stop_after:
                self.halt_on_next = True

        elapsed = self.time_elapsed(self.last_timestamp)
        t = time.strftime("%m-%d %H:%M:%S")
        if checkpoint is None:
            msg = "{m} ({t}) elapsed: {delta_t} _TIME_".format(
                m=message, t=t, delta_t=elapsed
            )
        else:
            msg = "{m} ({t}) ({status} {stage}) elapsed: {delta_t} _TIME_".format(
                m=message,
                t=t,
                status="finished" if finished else "starting",
                stage=checkpoint,
                delta_t=elapsed,
            )
        if re.match("^###", message):
            msg = "\n{}\n".format(msg)
        self.info(msg)
        self.last_timestamp = time.time()

    @staticmethod
    def time_elapsed(time_since):
        """
        Returns the number of seconds that have elapsed since the time_since parameter.

        :param float time_since: Time as a float given by time.time().
        """
        return round(time.time() - time_since, 0)

    def _report_profile(
        self, command, lock_name, elapsed_time, memory, pid, args_hash, proc_count
    ):
        """
        Writes a string to self.pipeline_profile_file.
        """
        rel_lock_name = (
            lock_name
            if lock_name is None
            else os.path.relpath(lock_name, self.outfolder)
        )
        message_raw = (
            str(pid)
            + "\t"
            + str(args_hash)
            + "\t"
            + str(proc_count)
            + "\t"
            + str(datetime.timedelta(seconds=round(elapsed_time, 2)))
            + "\t "
            + str(round(memory, 4))
            + "\t"
            + str(command)
            + "\t"
            + str(rel_lock_name)
        )
        with open(self.pipeline_profile_file, "a") as myfile:
            myfile.write(message_raw + "\n")

    def report_result(
        self, key, value, nolog=False, result_formatter=None, force_overwrite=True
    ):
        """
        Writes a key:value pair to self.pipeline_stats_file.

        :param str key: name (key) of the stat
        :param dict value: value of the stat to report.
        :param bool nolog: Turn on this flag to NOT print this result in the
            logfile. Use sparingly in case you will be printing the result in a
            different format.
        :param str result_formatter: function for formatting via pipestat backend
        :param bool force_overwrite: overwrite results if they already exist?
        :return str reported_result: the reported result is returned as a list of formatted strings.

        """
        # keep the value in memory:
        self.stats_dict[key] = value

        rf = result_formatter or self.pipestat_result_formatter

        reported_result = self.pipestat.report(
            values={key: value},
            record_identifier=self.pipestat_record_identifier,
            result_formatter=rf,
            force_overwrite=force_overwrite,
        )

        if not nolog:
            if isinstance(
                reported_result, bool
            ):  # Pipestat can return False if results are NOT reported.
                self.info("Result successfully reported? " + str(reported_result))
            else:
                for r in reported_result:
                    self.info(r)

        return reported_result

    def report_object(
        self,
        key,
        filename,
        anchor_text=None,
        anchor_image=None,
        annotation=None,
        nolog=False,
        result_formatter=None,
        force_overwrite=True,
    ):
        """
        Writes a key:value pair to self.pipeline_stats_file. Note: this function
            will be deprecated. Using report_result is recommended.

        :param str key: name (key) of the object
        :param str filename: relative path to the file (relative to parent
            output dir)
        :param str anchor_text: text used as the link anchor test or caption to
            refer to the object. If not provided, defaults to the key.
        :param str anchor_image: a path to an HTML-displayable image thumbnail
            (so, .png or .jpg, for example). If a path, the path should be
            relative to the parent output dir.
        :param str annotation: By default, the figures will be annotated with
            the pipeline name, so you can tell which pipeline records which
            figures. If you want, you can change this.
        :param bool nolog: Turn on this flag to NOT print this result in the
            logfile. Use sparingly in case you will be printing the result in a
            different format.
        :param str result_formatter: function for formatting via pipestat backend
        :param bool force_overwrite: overwrite results if they already exist?
        :return str reported_result: the reported result is returned as a list of formatted strings.
        """
        warnings.warn(
            "This function may be removed in future release. "
            "The recommended way to report pipeline results is using PipelineManager.pipestat.report().",
            category=DeprecationWarning,
        )
        rf = result_formatter or self.pipestat_result_formatter
        # Default annotation is current pipeline name.
        annotation = str(annotation or self.name)
        # In case the value is passed with trailing whitespace.
        filename = str(filename).strip()
        if anchor_text:
            anchor_text = str(anchor_text).strip()
        else:
            anchor_text = str(key).strip()
        # better to use a relative path in this file
        # convert any absolute paths into relative paths

        values = {
            "path": filename,
            "thumbnail_path": anchor_image,
            "title": anchor_text,
            "annotation": annotation,
        }
        val = {key: values}

        reported_result = self.pipestat.report(
            values=val,
            record_identifier=self.pipestat_record_identifier,
            result_formatter=rf,
            force_overwrite=force_overwrite,
        )

        if not nolog:
            if isinstance(
                reported_result, bool
            ):  # Pipestat can return False if results are NOT reported.
                self.info("Result successfully reported? " + str(reported_result))
            else:
                for r in reported_result:
                    self.info(r)

    def _safe_write_to_file(self, file, message):
        """
        Writes a string to a file safely (with file locks).
        """
        warnings.warn(
            "This function may be removed in future release. "
            "The recommended way to report pipeline results is using PipelineManager.pipestat.report().",
            category=DeprecationWarning,
        )
        target = file
        lock_name = make_lock_name(target, self.outfolder)
        lock_file = self._make_lock_path(lock_name)

        while True:
            if os.path.isfile(lock_file):
                self._wait_for_lock(lock_file)
            else:
                try:
                    self.locks.append(lock_file)
                    self._create_file_racefree(lock_file)
                except OSError as e:
                    if e.errno == errno.EEXIST:
                        self.warning("Lock file created after test! Looping again.")
                        continue  # Go back to start

                # Proceed with file writing
                with open(file, "a") as myfile:
                    myfile.write(message + "\n")

                os.remove(lock_file)
                self.locks.remove(lock_file)

                # If you make it to the end of the while loop, you're done
                break

    def _report_command(self, cmd, procs=None):
        """
        Writes a command to both stdout and to the commands log file
        (self.pipeline_commands_file).

        :param str cmd: command to report
        :param str | list[str] procs: process numbers for processes in the command
        """
        if isinstance(procs, list):
            procs = ",".join(map(str, procs))
        if procs:
            line = "\n> `{cmd}` ({procs})".format(cmd=str(cmd), procs=procs)
        else:
            line = "\n> `{cmd}`".format(cmd=str(cmd))

        # Print line to stdout
        self.info(line)

        # And also add to commands file

        cmd_line = "{cmd}\n".format(cmd=str(cmd))
        with open(self.pipeline_commands_file, "a") as myfile:
            myfile.write(cmd_line)

    ###################################
    # Filepath functions
    ###################################

    @staticmethod
    def _create_file(file):
        """
        Creates a file, but will not fail if the file already exists.
        This is vulnerable to race conditions; use this for cases where it
        doesn't matter if this process is the one that created the file.

        :param str file: File to create.
        """
        with open(file, "w") as fout:
            fout.write("")

    @staticmethod
    def _create_file_racefree(file):
        """
        Creates a file, but fails if the file already exists.

        This function will thus only succeed if this process actually creates
        the file; if the file already exists, it will cause an OSError,
        solving race conditions.

        :param str file: File to create.
        """
        write_lock_flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY
        fd = os.open(file, write_lock_flags)
        os.close(fd)

    @staticmethod
    def _ensure_lock_prefix(lock_name_base):
        """Ensure that an alleged lock file is correctly prefixed."""
        return (
            lock_name_base
            if lock_name_base.startswith(LOCK_PREFIX)
            else LOCK_PREFIX + lock_name_base
        )

    def _make_lock_path(self, lock_name_base):
        """
        Create path to lock file with given name as base.

        :param str lock_name_base: Lock file name, designed to not be prefixed
            with the lock file designation, but that's permitted.
        :return str: Path to the lock file.
        """

        # For lock prefix validation, separate file name from other path
        # components, as we care about the name prefix not path prefix.
        base, name = os.path.split(lock_name_base)

        lock_name = self._ensure_lock_prefix(name)
        if base:
            lock_name = os.path.join(base, lock_name)
        return pipeline_filepath(self, filename=lock_name)

    def _recoverfile_from_lockfile(self, lockfile):
        """
        Create path to recovery file with given name as base.

        :param str lockfile: Name of file on which to base this path,
            perhaps already prefixed with the designation of a lock file.
        :return str: Path to recovery file.
        """
        # Require that the lock file path be absolute, or at least relative
        # and starting with the pipeline output folder.
        if not (os.path.isabs(lockfile) or lockfile.startswith(self.outfolder)):
            lockfile = self._make_lock_path(lockfile)
        return lockfile.replace(LOCK_PREFIX, "recover." + LOCK_PREFIX)

    @staticmethod
    def make_sure_path_exists(path):
        """
        Creates all directories in a path if it does not exist.

        :param str path: Path to create.
        :raises Exception: if the path creation attempt hits an error with
            a code indicating a cause other than pre-existence.
        """
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    ###################################
    # Pipeline stats calculation helpers
    ###################################

    def _refresh_stats(self):
        """
        Loads up the stats yaml created for this pipeline run and reads
        those stats into memory
        """

        if os.path.isfile(self.pipeline_stats_file):
            data = load_yaml(filepath=self.pipeline_stats_file)

            for key, value in data[self._pipestat_manager.pipeline_name][
                self._pipestat_manager.pipeline_type
            ][self._pipestat_manager.record_identifier].items():
                self.stats_dict[key] = value

    def get_stat(self, key):
        """
        Returns a stat that was previously reported. This is necessary for reporting new stats that are
        derived from two stats, one of which may have been reported by an earlier run. For example,
        if you first use report_result to report (number of trimmed reads), and then in a later stage
        want to report alignment rate, then this second stat (alignment rate) will require knowing the
        first stat (number of trimmed reads); however, that may not have been calculated in the current
        pipeline run, so we must retrieve it from the stats.yaml output file. This command will retrieve
        such previously reported stats if they were not already calculated in the current pipeline run.
        :param key: key of stat to retrieve
        """

        try:
            return self.stats_dict[key]
        except KeyError:
            self._refresh_stats()
            try:
                return self.stats_dict[key]
            except KeyError:
                self.warning("Missing stat '{}'".format(key))
                return None

    ###################################
    # Pipeline termination functions
    ###################################

    # TODO: support checkpointing via stage, basic name, function, or filename (prohibit filepath)

    def _checkpoint(self, stage):
        """
        Decide whether to stop processing of a pipeline. This is the hook

        A pipeline can report various "checkpoints" as sort of status markers
        that designate the logical processing phase that's just been completed.
        The initiation of a pipeline can preordain one of those as a "stopping
        point" that when reached, should stop the pipeline's execution.

        :param pypiper.Stage | str stage: Pipeline processing stage/phase just completed.
        :return bool: Whether a checkpoint was created (i.e., whether it didn't
            already exist)
        :raise ValueError: If the stage is specified as an absolute filepath,
            and that path indicates a location that's not immediately within
            the main output folder, raise a ValueError.
        """

        # For null stage, short-circuit and indicate no file write.
        # This handles case in which we're timestamping prospectively and
        # previously weren't in a stage.
        if stage is None:
            return False

        try:
            is_checkpoint = stage.checkpoint
        except AttributeError:
            # Maybe we have a raw function, not a stage.
            if hasattr(stage, "__call__"):
                stage = stage.__name__
            else:
                # Maybe we have a stage name not a Stage.
                # In that case, we can proceed as-is, with downstream
                # processing handling Stage vs. stage name disambiguation.
                # Here, though, warn about inputs that appear filename/path-like.
                # We can't rely on raw text being a filepath or filename,
                # because that would ruin the ability to pass stage name rather
                # than actual stage. We can issue a warning message based on the
                # improbability of a stage name containing the '.' that would
                # be expected to characterize the extension of a file name/path.
                base, ext = os.path.splitext(stage)
                if ext and "." not in base:
                    self.warning(
                        "WARNING: '{}' looks like it may be the name or path of "
                        "a file; for such a checkpoint, use touch_checkpoint.".format(
                            stage
                        )
                    )
        else:
            if not is_checkpoint:
                self.warning("Not a checkpoint: {}".format(stage))
                return False
            stage = stage.name

        self.info("Checkpointing: '{}'".format(stage))
        if os.path.isabs(stage):
            check_fpath = stage
        else:
            check_fpath = checkpoint_filepath(stage, pm=self)
        return self._touch_checkpoint(check_fpath)

    def _touch_checkpoint(self, check_file):
        """
        Alternative way for a pipeline to designate a checkpoint.

        :param str check_file: Name or path of file to use as checkpoint.
        :return bool: Whether a file was written (equivalent to whether the
            checkpoint file already existed).
        :raise ValueError: Raise a ValueError if the argument provided as the
            checkpoint file is an absolute path and that doesn't correspond
            to a location within the main output folder.
        """
        if os.path.isabs(check_file):
            folder, _ = os.path.split(check_file)
            # For raw string comparison, ensure that each path
            # bears the final path separator.
            other_folder = os.path.join(folder, "")
            this_folder = os.path.join(self.outfolder, "")
            if other_folder != this_folder:
                errmsg = (
                    "Path provided as checkpoint file isn't in pipeline "
                    "output folder. '{}' is not in '{}'".format(
                        check_file, self.outfolder
                    )
                )
                raise ValueError(errmsg)
            fpath = check_file
        else:
            fpath = pipeline_filepath(self, filename=check_file)

        # Create/update timestamp for checkpoint, but base return value on
        # whether the action was a simple update or a novel creation.
        already_exists = os.path.isfile(fpath)
        open(fpath, "w").close()
        action = "Updated" if already_exists else "Created"
        self.info("{} checkpoint file: '{}'".format(action, fpath))

        return already_exists

    def complete(self):
        """Stop a completely finished pipeline."""
        self.stop_pipeline(status=COMPLETE_FLAG)

    def fail_pipeline(self, exc: Exception, dynamic_recover: bool = False):
        """
        If the pipeline does not complete, this function will stop the pipeline gracefully.
        It sets the status flag to failed and skips the normal success completion procedure.

        :param Exception exc: Exception to raise.
        :param bool dynamic_recover: Whether to recover e.g. for job termination.
        """
        # Take care of any active running subprocess
        sys.stdout.flush()
        self._terminate_running_subprocesses()

        if dynamic_recover:
            # job was terminated, not failed due to a bad process.
            # flag this run as recoverable.
            if len(self.locks) < 1:
                # If there is no process locked, then recovery will be automatic.
                self.info("No locked process. Dynamic recovery will be automatic.")
            # make a copy of self.locks to iterate over since we'll be clearing them as we go
            # set a recovery flag for each lock.
            for lock_file in self.locks[:]:
                recover_file = self._recoverfile_from_lockfile(lock_file)
                self.info("Setting dynamic recover file: {}".format(recover_file))
                self._create_file(recover_file)
                self.locks.remove(lock_file)

        # Produce cleanup script
        self._cleanup(dry_run=True)

        # Finally, set the status to failed and close out with a timestamp
        if not self._failed:  # and not self._completed:
            self.timestamp("### Pipeline failed at: ")
            total_time = datetime.timedelta(seconds=self.time_elapsed(self.starttime))
            self.info("Total time: " + str(total_time))
            self.info("Failure reason: " + str(exc))
            self.pipestat.set_status(
                record_identifier=self._pipestat_manager.record_identifier,
                status_identifier="failed",
            )

        if isinstance(exc, str):
            exc = RuntimeError(exc)

        raise exc

    def halt(self, checkpoint=None, finished=False, raise_error=True):
        """
        Stop the pipeline before completion point.

        :param str checkpoint: Name of stage just reached or just completed.
        :param bool finished: Whether the indicated stage was just finished
            (True), or just reached (False)
        :param bool raise_error: Whether to raise an exception to truly
            halt execution.
        """
        self.stop_pipeline(PAUSE_FLAG)
        self._active = False
        if raise_error:
            raise PipelineHalt(checkpoint, finished)

    def get_elapsed_time(self):
        """
        Parse the pipeline profile file, collect the unique and last duplicated
        runtimes and sum them up. In case the profile is not found, an estimate
        is calculated (which is correct only in case the pipeline was not rerun)

        :return int: sum of runtimes in seconds
        """
        if os.path.isfile(self.pipeline_profile_file):
            df = _pd.read_csv(
                self.pipeline_profile_file,
                sep="\t",
                comment="#",
                names=PROFILE_COLNAMES,
            )
            try:
                df["runtime"] = _pd.to_timedelta(df["runtime"])
            except ValueError:
                # return runtime estimate
                # this happens if old profile style is mixed with the new one
                # and the columns do not match
                return self.time_elapsed(self.starttime)
            unique_df = df[~df.duplicated("cid", keep="last").values]
            return sum(unique_df["runtime"].apply(lambda x: x.total_seconds()))
        return self.time_elapsed(self.starttime)

    def stop_pipeline(self, status=COMPLETE_FLAG):
        """
        Terminate the pipeline.

        This is the "healthy" pipeline completion function.
        The normal pipeline completion function, to be run by the pipeline
        at the end of the script. It sets status flag to completed and records
        some time and memory statistics to the log file.
        """
        # self._set_status_flag(status)
        self.pipestat.set_status(
            record_identifier=self._pipestat_manager.record_identifier,
            status_identifier=status,
        )
        self._cleanup()
        elapsed_time_this_run = str(
            datetime.timedelta(seconds=self.time_elapsed(self.starttime))
        )
        self.report_result("Time", elapsed_time_this_run, nolog=True)
        self.report_result("Success", time.strftime("%m-%d-%H:%M:%S"), nolog=True)

        self.info("\n### Pipeline completed. Epilogue")
        # print("* " + "Total elapsed time".rjust(20) + ":  "
        #       + str(datetime.timedelta(seconds=self.time_elapsed(self.starttime))))
        self.info(
            "* " + "Elapsed time (this run)".rjust(30) + ":  " + elapsed_time_this_run
        )
        self.info(
            "* "
            + "Total elapsed time (all runs)".rjust(30)
            + ":  "
            + str(datetime.timedelta(seconds=round(self.get_elapsed_time())))
        )
        self.info(
            "* "
            + "Peak memory (this run)".rjust(30)
            + ":  "
            + str(round(self.peak_memory, 4))
            + " GB"
        )
        # self.info("* " + "Total peak memory (all runs)".rjust(30) + ":  " +
        #     str(round(self.peak_memory, 4)) + " GB")
        if self.halted:
            return

        t = time.strftime("%Y-%m-%d %H:%M:%S")
        self.info("* " + "Pipeline completed time".rjust(30) + ": " + t)

    def _signal_term_handler(self, signal, frame):
        """
        TERM signal handler function: this function is run if the process
        receives a termination signal (TERM). This may be invoked, for example,
        by SLURM if the job exceeds its memory or time limits. It will simply
        record a message in the log file, stating that the process was
        terminated, and then gracefully fail the pipeline. This is necessary to
        1. set the status flag and 2. provide a meaningful error message in the
        tee'd output; if you do not handle this, then the tee process will be
        terminated before the TERM error message, leading to a confusing log
        file.
        """
        signal_type = "SIGTERM"
        self._generic_signal_handler(signal_type)

    def _generic_signal_handler(self, signal_type):
        """
        Function for handling both SIGTERM and SIGINT
        """
        self.info("</pre>")
        message = "Got " + signal_type + ". Failing gracefully..."
        self.timestamp(message)
        self.fail_pipeline(KeyboardInterrupt(signal_type), dynamic_recover=True)
        sys.exit(1)

        # I used to write to the logfile directly, because the interrupts
        # would first destroy the tee subprocess, so anything printed
        # after an interrupt would only go to screen and not get put in
        # the log, which is bad for cluster processing. But I later
        # figured out how to sequester the kill signals so they didn't get
        # passed directly to the tee subprocess, so I could handle that on
        # my own; hence, now I believe I no longer need to do this. I'm
        # leaving this code here as a relic in case something comes up.
        # with open(self.pipeline_log_file, "a") as myfile:
        #   myfile.write(message + "\n")

    def _signal_int_handler(self, signal, frame):
        """
        For catching interrupt (Ctrl +C) signals. Fails gracefully.
        """
        signal_type = "SIGINT"
        self._generic_signal_handler(signal_type)

    def _exit_handler(self):
        """
        This function I register with atexit to run whenever the script is completing.
        A catch-all for uncaught exceptions, setting status flag file to failed.
        """

        # TODO: consider handling sys.stderr/sys.stdout exceptions related to
        # TODO (cont.): order of interpreter vs. subprocess shutdown signal receipt.
        # TODO (cont.): see https://bugs.python.org/issue11380

        # Make the cleanup file executable if it exists
        if os.path.isfile(self.cleanup_file):
            # Make the cleanup file self destruct.
            with open(self.cleanup_file, "a") as myfile:
                myfile.write("rm " + self.cleanup_file + "\n")
            os.chmod(self.cleanup_file, 0o755)

        # If the pipeline hasn't completed successfully, or already been marked
        # as failed, then mark it as failed now.

        if not self._has_exit_status:
            self.debug("Pipeline status: {}".format(self.status))
            self.fail_pipeline(Exception("Pipeline failure. See details above."))

        if self.tee:
            self.tee.kill()

    def _terminate_running_subprocesses(self):
        # make a copy of the list to iterate over since we'll be removing items
        for pid in self.running_procs.copy():
            proc_dict = self.running_procs[pid]

            # Close the preformat tag that we opened when the process was spawned.
            # record profile of any running processes before killing
            elapsed_time = time.time() - self.running_procs[pid]["start_time"]
            process_peak_mem = (
                self._memory_usage(pid, container=proc_dict["container"]) / 1e6
            )
            self._report_profile(
                self.running_procs[pid]["proc_name"],
                None,
                elapsed_time,
                process_peak_mem,
                pid,
                self.running_procs[pid]["args_hash"],
                self.running_procs[pid]["local_proc_id"],
            )
            self._kill_child_process(pid, proc_dict["proc_name"])
            del self.running_procs[pid]

    def _kill_child_process(self, child_pid, proc_name=None):
        """
        Pypiper spawns subprocesses. We need to kill them to exit gracefully,
        in the event of a pipeline termination or interrupt signal.
        By default, child processes are not automatically killed when python
        terminates, so Pypiper must clean these up manually.
        Given a process ID, this function just kills it.

        :param int child_pid: Child process id.
        """

        # When we kill process, it turns into a zombie, and we have to reap it.
        # So we can't just kill it and then let it go; we call wait

        def pskill(proc_pid, sig=signal.SIGINT):
            parent_process = psutil.Process(proc_pid)
            for child_proc in parent_process.children(recursive=True):
                child_proc.send_signal(sig)
            parent_process.send_signal(sig)

        if child_pid is None:
            return

        if proc_name:
            proc_string = " ({proc_name})".format(proc_name=proc_name)

        # First a gentle kill
        sys.stdout.flush()
        still_running = self._attend_process(psutil.Process(child_pid), 0)
        sleeptime = 0.25
        time_waiting = 0

        while still_running and time_waiting < 3:
            try:
                if time_waiting > 2:
                    pskill(child_pid, signal.SIGKILL)
                    # print("pskill("+str(child_pid)+", signal.SIGKILL)")
                elif time_waiting > 1:
                    pskill(child_pid, signal.SIGTERM)
                    # print("pskill("+str(child_pid)+", signal.SIGTERM)")
                else:
                    pskill(child_pid, signal.SIGINT)
                    # print("pskill("+str(child_pid)+", signal.SIGINT)")

            except OSError:
                # This would happen if the child process ended between the check
                # and the next kill step
                still_running = False
                time_waiting = time_waiting + sleeptime

            # Now see if it's still running
            time_waiting = time_waiting + sleeptime
            if not self._attend_process(psutil.Process(child_pid), sleeptime):
                still_running = False

        if still_running:
            # still running!?
            self.warning(
                "Child process {child_pid}{proc_string} never responded"
                "I just can't take it anymore. I don't know what to do...".format(
                    child_pid=child_pid, proc_string=proc_string
                )
            )
        else:
            if time_waiting > 0:
                note = "terminated after {time} sec".format(time=int(time_waiting))
            else:
                note = "was already terminated"

            msg = "Child process {child_pid}{proc_string} {note}.".format(
                child_pid=child_pid, proc_string=proc_string, note=note
            )
            self.info(msg)

    @staticmethod
    def _atexit_register(*args):
        """Convenience alias to register exit functions without having to import atexit in the pipeline."""
        atexit.register(*args)

    def get_container(self, image, mounts):
        # image is something like "nsheff/refgenie"
        if type(mounts) == str:
            mounts = [mounts]
        cmd = "docker run -itd"
        for mnt in mounts:
            absmnt = os.path.abspath(mnt)
            cmd += " -v " + absmnt + ":" + absmnt
        cmd += " -v {cwd}:{cwd} --workdir={cwd}".format(cwd=os.getcwd())
        cmd += " --user={uid}".format(uid=os.getuid())
        cmd += " --volume=/etc/group:/etc/group:ro"
        cmd += " --volume=/etc/passwd:/etc/passwd:ro"
        cmd += " --volume=/etc/shadow:/etc/shadow:ro"
        cmd += " --volume=/etc/sudoers.d:/etc/sudoers.d:ro"
        cmd += " --volume=/tmp/.X11-unix:/tmp/.X11-unix:rw"
        cmd += " " + image
        container = self.checkprint(cmd).rstrip()
        self.container = container
        self.info("Using docker container: " + container)
        self._atexit_register(self.remove_container, container)

    def remove_container(self, container):
        self.info("Removing docker container. . .")
        cmd = "docker rm -f " + container
        self.callprint(cmd)

    def clean_add(self, regex, conditional=False, manual=False):
        """
        Add files (or regexs) to a cleanup list, to delete when this pipeline completes successfully.
        When making a call with run that produces intermediate files that should be
        deleted after the pipeline completes, you flag these files for deletion with this command.
        Files added with clean_add will only be deleted upon success of the pipeline.

        :param str regex:  A unix-style regular expression that matches files to delete
            (can also be a file name).
        :param bool conditional: True means the files will only be deleted if no other
            pipelines are currently running; otherwise they are added to a manual cleanup script
            called {pipeline_name}_cleanup.sh
        :param bool manual: True means the files will just be added to a manual cleanup script.
        """

        # TODO: print this message (and several below) in debug
        # print("Adding regex to cleanup: {}".format(regex))
        if self.dirty:
            # Override the user-provided option and force manual cleanup.
            manual = True

        if not self.clean_initialized:
            # Make cleanup files relative to the cleanup script in case the result folder moves.
            with open(self.cleanup_file, "a") as myfile:
                clean_init = 'DIR="$(cd -P -- "$(dirname -- "$0")" && pwd -P)"'
                myfile.write(clean_init + "\n")
                myfile.write("cd ${DIR}\n")
                self.clean_initialized = True

        if manual:
            filenames = glob.glob(regex)
            if not filenames:
                self.info("No files match cleanup pattern: {}".format(regex))
            for filename in filenames:
                try:
                    with open(self.cleanup_file, "a") as myfile:
                        if os.path.isabs(filename):
                            relative_filename = os.path.relpath(
                                filename, self.outfolder
                            )
                            absolute_filename = filename
                        else:
                            relative_filename = os.path.relpath(
                                filename, self.outfolder
                            )
                            absolute_filename = os.path.abspath(
                                os.path.join(self.outfolder, relative_filename)
                            )
                        if os.path.isfile(absolute_filename):
                            # print("Adding file to cleanup: {}".format(filename))
                            myfile.write("rm " + relative_filename + "\n")
                        elif os.path.isdir(absolute_filename):
                            # print("Adding directory to cleanup: {}".format(filename))
                            # first, add all filenames in the directory
                            myfile.write("rm " + relative_filename + "/*\n")
                            # and the directory itself
                            myfile.write("rmdir " + relative_filename + "\n")
                        else:
                            self.info(
                                "File not added to cleanup: {}".format(
                                    relative_filename
                                )
                            )
                except Exception as e:
                    self.error(
                        "Error in clean_add on path {}: {}".format(filename, str(e))
                    )
        elif conditional:
            self.cleanup_list_conditional.append(regex)
        else:
            self.cleanup_list.append(regex)
            # TODO: what's the "absolute" list?
            # Remove it from the conditional list if added to the absolute list
            while regex in self.cleanup_list_conditional:
                self.cleanup_list_conditional.remove(regex)

    def _cleanup(self, dry_run=False):
        """
        Cleans up (removes) intermediate files.

        You can register intermediate files, which will be deleted automatically
        when the pipeline completes. This function deletes them,
        either absolutely or conditionally. It is run automatically when the
        pipeline succeeds, so you shouldn't need to call it from a pipeline.

        :param bool dry_run: Set to True if you want to build a cleanup
            script, but not actually delete the files.
        """

        n_to_clean = len(self.cleanup_list)
        n_to_clean_cond = len(self.cleanup_list_conditional)

        if n_to_clean + n_to_clean_cond > 0:
            self.info(
                "Starting cleanup: {} files; {} conditional files for cleanup".format(
                    n_to_clean, n_to_clean_cond
                )
            )
        else:
            self.debug("No files to clean.")

        if dry_run:
            # Move all unconditional cleans into the conditional list
            if n_to_clean > 0:
                combined_list = self.cleanup_list_conditional + self.cleanup_list
                self.cleanup_list_conditional = combined_list
                self.cleanup_list = []

        if n_to_clean > 0:
            self.info("\nCleaning up flagged intermediate files. . .")
            for expr in self.cleanup_list:
                self.debug("Removing glob: " + expr)
                try:
                    # Expand regular expression
                    files = glob.glob(expr)
                    # Remove entry from cleanup list
                    while files in self.cleanup_list:
                        self.cleanup_list.remove(files)
                    # and delete the files
                    for file in files:
                        if os.path.isfile(file):
                            self.debug("`rm {}`".format(file))
                            os.remove(os.path.join(file))
                        elif os.path.isdir(file):
                            self.debug("`rmdir {}`".format(file))
                            os.rmdir(os.path.join(file))
                except:
                    pass

        if n_to_clean_cond > 0:
            run_flag = flag_name(RUN_FLAG)
            flag_files = [
                fn
                for fn in glob.glob(self.outfolder + flag_name("*"))
                if COMPLETE_FLAG not in os.path.basename(fn)
                and not "{}_{}_{}".format(
                    self._pipestat_manager.pipeline_name,
                    self.pipestat_record_identifier,
                    run_flag,
                )
                == os.path.basename(fn)
            ]
            if len(flag_files) == 0 and not dry_run:
                self.info("\nCleaning up conditional list. . .")
                for expr in self.cleanup_list_conditional:
                    self.debug("\nRemoving glob: " + expr)
                    try:
                        files = glob.glob(expr)
                        while files in self.cleanup_list_conditional:
                            self.cleanup_list_conditional.remove(files)
                        for file in files:
                            if os.path.isfile(file):
                                self.debug("`rm {}`".format(file))
                                os.remove(os.path.join(file))
                            elif os.path.isdir(file):
                                self.debug("`rmdir {}`".format(file))
                                os.rmdir(os.path.join(file))
                    except:
                        pass
            else:
                self.info(
                    "\nConditional flag found: "
                    + str([os.path.basename(i) for i in flag_files])
                )
                self.info(
                    "\nThese conditional files were left in place:\n\n- "
                    + "\n- ".join(self.cleanup_list_conditional)
                )
                # Produce a cleanup script.
                no_cleanup_script = []
                for cleandir in self.cleanup_list_conditional:
                    try:
                        items_to_clean = glob.glob(cleandir)
                        for clean_item in items_to_clean:
                            with open(self.cleanup_file, "a") as clean_script:
                                if os.path.isfile(clean_item):
                                    clean_script.write("rm " + clean_item + "\n")
                                elif os.path.isdir(clean_item):
                                    clean_script.write("rmdir " + clean_item + "\n")
                    except Exception as e:
                        no_cleanup_script.append(cleandir)
                if no_cleanup_script:
                    self.warning(
                        "\n\nCould not produce cleanup script for item(s):\n\n- "
                        + "\n- ".join(no_cleanup_script)
                    )

    def _memory_usage(self, pid="self", category="hwm", container=None):
        """
        Memory usage of the process in kilobytes.

        :param str pid: Process ID of process to check
        :param str category: Memory type to check. 'hwm' for high water mark.
        """
        if container:
            # TODO: Put some debug output here with switch to Logger
            # since this is relatively untested.
            cmd = "docker stats " + container + " --format '{{.MemUsage}}' --no-stream"
            mem_use_str = subprocess.check_output(cmd, shell=True).decode()

            mem_num = re.findall(r"[\d\.]+", mem_use_str.split("/")[0])[0]
            mem_scale = re.findall(r"[A-Za-z]+", mem_use_str.split("/")[0])[0]

            mem_num = float(mem_num)
            if mem_scale == "GiB":
                return mem_num * 1e6
            elif mem_scale == "MiB":
                return mem_num * 1e3
            elif mem_scale == "KiB":
                return mem_num
            else:
                # What type is this?
                return 0

        # Thanks Martin Geisler:
        status = None
        result = {"peak": 0, "rss": 0, "hwm": 0}

        try:
            # This will only work on systems with a /proc file system
            # (like Linux).
            # status = open('/proc/self/status')
            proc_spot = "/proc/%s/status" % pid
            status = open(proc_spot)
            for line in status:
                parts = line.split()
                key = parts[0][2:-1].lower()
                if key in result:
                    result[key] = int(parts[1])
        except:
            return 0

        finally:
            if status is not None:
                status.close()
        # print(result[category])
        return result[category]

    def _triage_error(self, e, nofail):
        """Print a message and decide what to do about an error."""
        if not nofail:
            self.fail_pipeline(e)
        elif self._failed:
            self.info(
                "This is a nofail process, but the pipeline was terminated for other reasons, so we fail."
            )
            raise e
        else:
            self.error(e)
            self.error(
                "ERROR: Subprocess returned nonzero result, but pipeline is continuing because nofail=True"
            )
        # TODO: return nonzero, or something. . .?
