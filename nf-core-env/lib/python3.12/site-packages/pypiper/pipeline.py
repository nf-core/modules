"""Pipeline base class"""

import abc
import glob
import os
import sys

from collections.abc import Iterable, Mapping

from .exceptions import (
    IllegalPipelineDefinitionError,
    IllegalPipelineExecutionError,
    UnknownPipelineStageError,
)
from .manager import PipelineManager
from .stage import Stage
from .utils import (
    checkpoint_filepath,
    flag_name,
    parse_stage_name,
    translate_stage_name,
)

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


__all__ = ["Pipeline", "UnknownPipelineStageError"]


class Pipeline(object):
    """
    Generic pipeline framework.

    Note that either a PipelineManager or an output folder path is required
    to create the Pipeline. If both are provided, the output folder path is
    ignored and that of the PipelineManager that's passed is used instead.

    :param str name: Name for the pipeline; arbitrary, just used for messaging;
        this is required if and only if a manager is not provided.
    :param pypiper.PipelineManager manager: The pipeline manager to use for
        this pipeline; this is required if and only if name or output folder is
        not provided.
    :param str outfolder: Path to main output folder for this pipeline
    :param argparse.Namespace args: command-line options and arguments for PipelineManager
    :param Mapping pl_mgr_kwargs: Additional keyword arguments for pipeline manager.
    :raise TypeError: Either pipeline manager or output folder path must be
        provided, and either pipeline manager or name must be provided.
    :raise ValueError: Name of pipeline manager cannot be
    :raise pypiper.IllegalPipelineDefinitionError: Definition of collection of
        stages must be non-empty.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(
        self, name=None, manager=None, outfolder=None, args=None, **pl_mgr_kwargs
    ):
        super(Pipeline, self).__init__()
        try:
            self.name = name or manager.name
        except AttributeError:
            raise TypeError(
                "If a pipeline manager isn't provided to create "
                "{}, a name is required.".format(Pipeline.__name__)
            )
        else:
            if not self.name:
                raise ValueError(
                    "Invalid name, possible inferred from pipeline manager: "
                    "{} ({})".format(self.name, type(self.name))
                )

        # Determine the PipelineManager.
        if manager:
            self.manager = manager
            if outfolder:
                print(
                    "Ignoring explicit output folder ({}) and using that of "
                    "pipeline manager ({})".format(outfolder, manager.outfolder)
                )
            if name and name != manager.name:
                print(
                    "Warning: name for pipeline ('{}') doesn't match that "
                    "of the given manager ('{}')".format(name, manager.name)
                )
        elif outfolder:
            # We're guaranteed by the upfront exception block around
            # name setting that we'll either have set the name for this
            # instance to a non-null if we reach this point, and thus we're
            # protected from passing a null name argument to the pipeline
            # manager's constructor.
            self.manager = PipelineManager(
                self.name, outfolder, args=args, **pl_mgr_kwargs
            )
        else:
            raise TypeError(
                "To create a {} instance, 'manager' or 'outfolder' "
                "is required".format(self.__class__.__name__)
            )

        # Require that checkpoints be overwritten.
        self.manager.overwrite_checkpoints = True

        # Translate stage names; do this here to complicate a hypothetical
        # attempt to override or otherwise redefine the way in which
        # stage names are handled, parsed, and translated.
        self._unordered = _is_unordered(self.stages())
        if self._unordered:
            print(
                "NOTICE: Unordered definition of stages for "
                "pipeline {}".format(self.name)
            )

        # Get to a sequence of pairs of key (possibly in need of translation)
        # and actual callable. Key is stage name and value is either stage
        # callable or an already-made stage object.
        stages = (
            self.stages().items()
            if isinstance(self.stages(), Mapping)
            else self.stages()
        )
        # Stage spec. parser handles callable validation.
        name_stage_pairs = [_parse_stage_spec(s) for s in stages]

        # Pipeline must have non-empty definition of stages.
        if not stages:
            raise IllegalPipelineDefinitionError("Empty stages")

        # Ensure that each pipeline stage is callable, and map names
        # between from external specification and internal representation.
        # We don't need to store the internal-to-external mapping, as each
        # stage will store its own name that is equivalent to the "external"
        # one, and we can use the checkpoint name derivation functions
        # to determine checkpoint name/path from stage/stage name.
        # We just use this internal-to-external mapping here, ephemerally,
        # to pretest whether there'd be a checkpoint name resolution collision.
        _internal_to_external = dict()
        self._external_to_internal = dict()
        self._stages = []

        for name, stage in name_stage_pairs:
            # Use external translator to further confound redefinition.
            internal_name = translate_stage_name(name)

            # Check that there's not a checkpoint name collision.
            if internal_name in _internal_to_external:
                already_mapped = _internal_to_external[internal_name]
                errmsg = (
                    "Duplicate stage name resolution (stage names are too "
                    "similar.) '{}' and '{}' both resolve to '{}'".format(
                        name, already_mapped, internal_name
                    )
                )
                raise IllegalPipelineDefinitionError(errmsg)

            # Store the stage name translations and the stage itself.
            self._external_to_internal[name] = internal_name
            _internal_to_external[internal_name] = name
            self._stages.append(stage)

        self.skipped, self.executed = None, None

    @property
    def outfolder(self):
        """
        Determine the path to the output folder for this pipeline instance.

        :return str: Path to output folder for this pipeline instance.
        """
        return self.manager.outfolder

    @abc.abstractmethod
    def stages(self):
        """
        Define the names of pipeline processing stages.

        :return Iterable[str]: Collection of pipeline stage names.
        """
        pass

    @property
    def stage_names(self):
        """
        Fetch the pipeline's stage names as specified by the pipeline
        class author (i.e., not necessarily those that are used for the
        checkpoint files)

        :return list[str]: Sequence of names of this pipeline's defined stages.
        """
        return [parse_stage_name(s) for s in self._stages]

    def checkpoint(self, stage, msg=""):
        """
        Touch checkpoint file for given stage and provide timestamp message.

        :param pypiper.Stage stage: Stage for which to mark checkpoint
        :param str msg: Message to embed in timestamp.
        :return bool: Whether a checkpoint file was written.
        """
        # Canonical usage model for Pipeline checkpointing through
        # implementations of this class is by automatically creating a
        # checkpoint when a conceptual unit or group of operations of a
        # pipeline completes, so fix the 'finished' parameter to the manager's
        # timestamp method to be True.
        return self.manager.timestamp(
            message=msg, checkpoint=stage.checkpoint_name, finished=True
        )

    def completed_stage(self, stage):
        """
        Determine whether the pipeline's completed the stage indicated.

        :param pypiper.Stage stage: Stage to check for completion status.
        :return bool: Whether this pipeline's completed the indicated stage.
        :raises UnknownStageException: If the stage name given is undefined
            for the pipeline, a ValueError arises.
        """
        check_path = checkpoint_filepath(stage, self.manager)
        return os.path.exists(check_path)

    def halt(self, **kwargs):
        """Halt the pipeline"""
        self.manager.halt(**kwargs)

    def list_flags(self, only_name=False):
        """
        Determine the flag files associated with this pipeline.

        :param bool only_name: Whether to return only flag file name(s) (True),
            or full flag file paths (False); default False (paths)
        :return list[str]: flag files associated with this pipeline.
        """
        paths = glob.glob(os.path.join(self.outfolder, flag_name("*")))
        if only_name:
            return [os.path.split(p)[1] for p in paths]
        else:
            return paths

    def run(self, start_point=None, stop_before=None, stop_after=None):
        """
        Run the pipeline, optionally specifying start and/or stop points.

        :param str start_point: Name of stage at which to begin execution.
        :param str stop_before: Name of stage at which to cease execution;
            exclusive, i.e. this stage is not run
        :param str stop_after: Name of stage at which to cease execution;
            inclusive, i.e. this stage is the last one run
        :raise IllegalPipelineExecutionError: If both inclusive (stop_after)
            and exclusive (stop_before) halting points are provided, or if that
            start stage is the same as or after the stop stage, raise an
            IllegalPipelineExecutionError.
        """

        # Start the run with a clean slate of Stage status/label tracking.
        self._reset()

        # TODO: validate starting point against checkpoint flags for
        # TODO (cont.): earlier stages if the pipeline defines its stages as a
        # TODO (cont.): sequence (i.e., probably prohibit start point with
        # TODO (cont): nonexistent earlier checkpoint flag(s).)

        if stop_before and stop_after:
            raise IllegalPipelineExecutionError(
                "Cannot specify both inclusive and exclusive stops."
            )

        if stop_before:
            stop = stop_before
            inclusive_stop = False
        elif stop_after:
            stop = stop_after
            inclusive_stop = True
        else:
            stop = None
            inclusive_stop = None

        # Ensure that a stage name--if specified--is supported.
        for s in [start_point, stop]:
            if s is None:
                continue
            name = parse_stage_name(s)
            if name not in self.stage_names:
                raise UnknownPipelineStageError(name, self)

        # Permit order-agnostic pipelines, but warn.
        if self._unordered and (start_point or stop_before or stop_after):
            print(
                "WARNING: Starting and stopping points are nonsense for "
                "pipeline with unordered stages."
            )

        # TODO: consider context manager based on start/stop points.

        # Determine where to start (but perhaps skip further based on
        # checkpoint completions.)
        start_index = self._start_index(start_point)
        stop_index = self._stop_index(stop, inclusive=inclusive_stop)
        assert stop_index <= len(self._stages)
        if start_index >= stop_index:
            raise IllegalPipelineExecutionError(
                "Cannot start pipeline at or after stopping point"
            )

        # TODO: consider storing just stage name rather than entire stage.
        # TODO (cont.): the bad case for whole-Stage is if associated data
        # TODO (cont.): (i.e., one or more args) are large.
        self.skipped.extend(self._stages[:start_index])

        # TODO: support both courses of action for non-continuous checkpoints.
        # TODO (cont.): That is, what if there's a stage with a checkpoint
        # TODO (cont.): file downstream of one without it? Naively, we'll
        # TODO (cont.): skip it, but we may want to re-run.
        skip_mode = True

        for stage in self._stages[start_index:stop_index]:
            # TODO: Note that there's no way to tell whether a non-checkpointed
            # TODO (cont.) Stage has been completed, and thus this seek
            # TODO (cont.) operation will find the first Stage, starting
            # TODO (cont.) the specified start point, either uncheckpointed or
            # TODO (cont.) for which the checkpoint file does not exist.
            # Look for checkpoint file.
            if skip_mode and self.completed_stage(stage):
                print("Skipping completed checkpoint stage: {}".format(stage))
                self.skipped.append(stage)
                continue

            # Once we've found where to being execution, ignore checkpoint
            # flags downstream if they exist since there may be dependence
            # between results from different stages.
            skip_mode = False

            print(f"Running stage: {getattr(stage, 'name', str(stage))}")

            try:
                stage.run()
            except Exception as e:
                self.manager._triage_error(e, nofail=stage.nofail)
            else:
                self.checkpoint(stage)
            self.executed.append(stage)

        # Add any unused stages to the collection of skips.
        self.skipped.extend(self._stages[stop_index:])

        # Where we stopped determines the shutdown mode.
        if stop_index == len(self._stages):
            self.wrapup()
        else:
            self.halt(raise_error=False)

    def wrapup(self):
        """Final mock stage to run after final one finishes."""
        self.manager.complete()

    def _reset(self):
        """Scrub decks with respect to Stage status/label tracking."""
        self.skipped, self.executed = [], []

    def _start_index(self, start=None):
        """Seek to the first stage to run."""
        if start is None:
            return 0
        start_stage = translate_stage_name(start)
        internal_names = [translate_stage_name(s.name) for s in self._stages]
        try:
            return internal_names.index(start_stage)
        except ValueError:
            raise UnknownPipelineStageError(start, self)

    def _stop_index(self, stop_point, inclusive):
        """
        Determine index of stage of stopping point for run().

        :param str | pypiper.Stage | function stop_point: Stopping point itself
            or name of it.
        :param bool inclusive: Whether the stopping point is to be regarded as
            inclusive (i.e., whether it's the final stage to run, or the one
            just beyond)
        :return int: Index into sequence of Pipeline's stages that indicates
            where to stop; critically, the value of the inclusive parameter
            here is used to contextualize this index such that it's always
            returned as an exclusive stopping index (i.e., execute up to the
            stage indexed by the value returned from this function.)
        """
        if not stop_point:
            # Null case, no stopping point
            return len(self._stages)
        stop_name = parse_stage_name(stop_point)
        try:
            stop_index = self.stage_names.index(stop_name)
        except ValueError:
            raise UnknownPipelineStageError(stop_name, self)
        return stop_index + 1 if inclusive else stop_index


def _is_unordered(collection):
    """
    Determine whether a collection appears to be unordered.

    This is a conservative implementation, allowing for the possibility that
    someone's implemented Mapping or Set, for example, and provided an
    __iter__ implementation that defines a consistent ordering of the
    collection's elements.

    :param object collection: Object to check as an unordered collection.
    :return bool: Whether the given object appears to be unordered
    :raises TypeError: If the given "collection" is non-iterable, it's
        illogical to investigate whether it's ordered.
    """
    if not isinstance(collection, Iterable):
        raise TypeError("Non-iterable alleged collection: {}".format(type(collection)))

    return isinstance(collection, set) or isinstance(collection, dict)


def _parse_stage_spec(stage_spec):
    """
    Handle alternate Stage specifications, returning name and Stage.

    Isolate this parsing logic from any iteration. TypeError as single
    exception type funnel also provides a more uniform way for callers to
    handle specification errors (e.g., skip a stage, warn, re-raise, etc.)

    :param (str, pypiper.Stage) | callable stage_spec: name and Stage
    :return (name, pypiper.Stage): Pair of name and Stage instance from parsing
        input specification
    :raise TypeError: if the specification of the stage is not a supported type
    """

    # The logic used here, a message to a user about how to specify Stage.
    req_msg = (
        "Stage specification must be either a {0} itself, a "
        "(<name>, {0}) pair, or a callable with a __name__ attribute "
        "(e.g., a non-anonymous function)".format(Stage.__name__)
    )

    # Simplest case is stage itself.
    if isinstance(stage_spec, Stage):
        return stage_spec.name, stage_spec

    # Handle alternate forms of specification.
    try:
        # Unpack pair of name and stage, requiring name first.
        name, stage = stage_spec
    except (TypeError, ValueError):
        # Normally, this sort of unpacking issue create a ValueError. Here,
        # though, we also need to catch TypeError since that's what arises
        # if an attempt is made to unpack a single function.
        # Attempt to parse stage_spec as a single named callable.
        try:
            name = stage_spec.__name__
        except AttributeError:
            raise TypeError(req_msg)
        else:
            # Control flow here indicates an anonymous function that was not
            # paired with a name. Prohibit that.
            if name == (lambda: None).__name__:
                raise TypeError(req_msg)
        stage = stage_spec

    # Ensure that the stage is callable.
    if not hasattr(stage, "__call__"):
        raise TypeError(req_msg)

    return name, Stage(stage, name=name)
