"""Project logging configuration.

Project-scope constants may reside here, but more importantly, some setup here
will provide a logging infrastructure for all of the project's modules.
Individual modules and classes may provide separate configuration on a more
local level, but this will at least provide a foundation.

"""

import logging
import os
import sys
import warnings
from ._version import __version__

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"

__all__ = [
    "add_logging_options",
    "logger_via_cli",
    "init_logger",
    "setup_logger",
    "AbsentOptionException",
    "LOGGING_CLI_OPTDATA",
]


BASIC_LOGGING_FORMAT = "%(message)s"
DEV_LOGGING_FMT = (
    "%(levelname).4s %(asctime)s | %(name)s:%(module)s:%(lineno)d > %(message)s "
)
FULL_DEV_LOGGING_FMT = (
    "%(levelname)s %(asctime)s | %(name)s:%(module)s:%(lineno)d > %(message)s "
)
DEFAULT_DATE_FMT = "%H:%M:%S"
PACKAGE_NAME = "logmuse"
STREAMS = {"OUT": sys.stdout, "ERR": sys.stderr}
DEFAULT_STREAM = STREAMS["ERR"]
LOGGING_LEVEL = "INFO"
LOGGING_LOCATIONS = (DEFAULT_STREAM,)
TRACE_LEVEL_VALUE = 5
TRACE_LEVEL_NAME = "TRACE"
CUSTOM_LEVELS = {TRACE_LEVEL_NAME: TRACE_LEVEL_VALUE}
SILENCE_LOGS_OPTNAME = "silent"
VERBOSITY_OPTNAME = "verbosity"
DEVMODE_OPTNAME = "logdev"
PARAM_BY_OPTNAME = {DEVMODE_OPTNAME: "devmode"}

# Translation of verbosity into logging level.
# Log message count monotonically increases in verbosity while it decreases
# in logging level, making verbosity a more intuitive specification mechanism.
_WARN_REPR = "WARN"
LEVEL_BY_VERBOSITY = ["CRITICAL", "ERROR", _WARN_REPR, "INFO", "DEBUG"]
_MIN_VERBOSITY = 1
_MAX_VERBOSITY = len(LEVEL_BY_VERBOSITY)
_VERBOSITY_CHOICES = (
    [str(x) for x in range(_MIN_VERBOSITY, len(LEVEL_BY_VERBOSITY) + 1)]
    + LEVEL_BY_VERBOSITY
    + ["WARNING"]
)

LOGGING_CLI_OPTDATA = {
    SILENCE_LOGS_OPTNAME: {
        "action": "store_true",
        "help": "Silence logging. Overrides {}.".format(VERBOSITY_OPTNAME),
    },
    VERBOSITY_OPTNAME: {
        "metavar": "V",
        "choices": _VERBOSITY_CHOICES,
        "help": "Set logging level ({}-{} or logging module level name)".format(
            _MIN_VERBOSITY, len(LEVEL_BY_VERBOSITY)
        ),
    },
    DEVMODE_OPTNAME: {
        "action": "store_true",
        "help": "Expand content of logging message format.",
    },
}


def add_logging_options(parser):
    """
    Augment a CLI argument parser with this package's logging options.

    :param argparse.ArgumentParser parser: CLI options and argument parser to
        augment with logging options.
    :return argparse.ArgumentParser: the input argument, supplemented with this
        package's logging options.
    """
    for optname, optdata in LOGGING_CLI_OPTDATA.items():
        parser.add_argument("--{}".format(optname), **optdata)
    return parser


def logger_via_cli(opts, strict=True, **kwargs):
    """
    Convenience function creating a logger.

    This module provides the ability to augment a CLI parser with
    logging-related options/arguments so that client applications do not need
    intimate knowledge of the implementation. This function completes that
    lack of burden, parsing values for the options supplied herein.

    :param argparse.Namespace opts: command-line options/arguments.
    :param bool strict: whether to raise an exception
    :return logging.Logger: configured logger instance.
    :raise pararead.logs.AbsentOptionException: if one of the expected options
        isn't available in the given Namespace, and the argument to the strict
        parameter is True. Such a case suggests that a client application
        didn't use this module to add the expected logging options to a parser.
    """
    # Within the key, translate the option name if needed. If it's not
    # present within the translations mapping, use the original optname.
    # Once translation's done (if needed), parse out the
    logs_cli_args = {}
    for optname in LOGGING_CLI_OPTDATA.keys():
        name = optname.lstrip("-")
        # Client must add the expected options, via the API or otherwise.
        try:
            optval = getattr(opts, name)
        except AttributeError:
            if strict:
                raise AbsentOptionException(optname)
            continue
        else:
            # Translate the option name if needed (i.e., for discordance
            # between the CLI version and the logger setup signature).
            logs_cli_args[PARAM_BY_OPTNAME.get(optname, name)] = optval
    logs_cli_args.update(kwargs)
    return init_logger(**logs_cli_args)


def init_logger(
    name="",
    level=None,
    stream=None,
    logfile=None,
    make_root=None,
    propagate=False,
    silent=False,
    devmode=False,
    verbosity=None,
    fmt=None,
    datefmt=DEFAULT_DATE_FMT,
    plain_format=False,
    style=None,
    use_full_names=False,
):
    """
    Establish and configure primary logger.

    This is intended to be called just once per "session", with a "session"
    defined as an invocation of the main workflow, a testing session, or an
    import of the primary abstractions, e.g. in an interactive iPython session.

    :param str name: name for the logger
    :param int | str level: minimal level of messages to listen for
    :param str stream: standard stream to use as log destination. The default
        behavior is to write logs to stdout, even if null is passed here. This
        is to allow a CLI argument as input to stream parameter, where it may be
        undesirable to require specification of a default value in the client
        application in order to prevent passing None if no CLI option value
        is given. To disable standard stream logging, set 'silent' to True
        or pass a path to a file to which to write logs, which gets priority
        over a standard stream as the destination for log messages.
    :param str | FileIO[str] logfile: path to filesystem location to use as
        logs destination. if provided, this mutes standard stream logging.
    :param bool make_root: whether to use returned logger as root logger. This
        means the name will be 'root' and that messages will not propagate.
    :param bool propagate: whether to allow messages from this logger to reach
        parent logger(s).
    :param bool silent: whether to silence logging; this is only guaranteed for
        messages from this logger and for those from loggers beneath this one
        in the runtime hierarchy without no separate handling. Propagation must
        also be turned off separately--if this is not the root logger--in
        order to ensure that messages are not handled and emitted from a
        potential parent to the logger built here.
    :param bool devmode: whether to log in development mode; possibly among
        other behavioral changes to logs handling, use a more information-rich
        message format template.
    :param int | str verbosity: alternate mode of expression for logging level
        that better accords with intuition about how to convey this. It's
        positively associated with message volume rather than negatively so, as
        logging level is. This takes precedence over 'level' if both are present.
    :param str fmt: message format/template.
    :param str datefmt: format/template for time component of a log record.
    :param bool plain_format: force use of plain message format, even if
        in development mode (debug level)
    :param str style: string indicating message formatting strategy; refer to
        https://docs.python.org/3/howto/logging-cookbook.html#use-of-alternative-formatting-styles;
        only valid in Python3.2+
    :param bool use_full_names: don't truncate level names
    :return logging.Logger: configured Logger instance
    :raise ValueError: if attempting to name explicitly non-root logger with
        a root name, or if both level and verbosity are specified
    """

    if make_root is True:
        if propagate:
            logging.warning("Propagation from root logger is nonsense")
        if name and name != "root":
            logging.warning(
                "Requested root logger with non-root name: " "{}".format(name)
            )
    else:
        name = name or PACKAGE_NAME
        if make_root is False and name == "root":
            raise ValueError(
                "Requested non-root logger with root name: {}".format(name)
            )

    # Enable named ultrafine logging for debugging.
    for level_name, level_value in CUSTOM_LEVELS.items():
        logging.addLevelName(level_value, level_name)

    # Establish the logger.
    logger = logging.getLogger(name)
    logger.handlers = []
    logger.propagate = propagate

    # Either short-circuit with a silent logger or parse and set level.
    if silent:
        logger.addHandler(logging.NullHandler())
        return logger

    # Determine the logger's listening level.
    if level is not None and verbosity is not None:
        raise ValueError(
            "Cannot specify both level and verbosity; got {} and "
            "{}, respectively".format(level, verbosity)
        )
    elif level is not None:
        # Handle int- or text-specific logging level.
        try:
            level = int(level)
        except ValueError:
            level = level.upper()
    else:
        level = _level_from_verbosity(verbosity or LOGGING_LEVEL)
    try:
        level = getattr(logging, level) if isinstance(level, str) else level
        logger.setLevel(level)
    except Exception:
        logging.error(
            "Can't set logging level to %s; instead using: '%s'",
            str(level),
            str(LOGGING_LEVEL),
        )
        level = LOGGING_LEVEL
        logger.setLevel(level)

    handlers = []

    if logfile:
        logfile_folder = os.path.dirname(logfile)
        if not os.path.exists(logfile_folder):
            os.makedirs(logfile_folder)

        # Create and add the handler, overwriting rather than appending.
        handlers.append(logging.FileHandler(logfile, mode="w"))
    if stream or not logfile:
        if not stream:
            stream = DEFAULT_STREAM
            stream_loc = stream
        elif stream in [sys.stderr, sys.stdout]:
            stream_loc = stream
        else:
            try:
                # Assume that we have a stream-indicative text argument.
                stream_loc = STREAMS[stream.upper()]
            except (AttributeError, KeyError):
                # Fall back on default stream since
                # arguments indicate that one should be activated.
                print(
                    "Invalid stream location: {}; using {}".format(
                        stream, DEFAULT_STREAM
                    )
                )
                stream_loc = DEFAULT_STREAM
        handlers.append(logging.StreamHandler(stream_loc))

    fine = level <= logging.DEBUG
    get_fmt = (
        (lambda _: fmt)
        if fmt
        else (
            lambda hdlr: (
                BASIC_LOGGING_FORMAT
                if plain_format
                or not (devmode or fine or isinstance(hdlr, logging.FileHandler))
                else (FULL_DEV_LOGGING_FMT if use_full_names else DEV_LOGGING_FMT)
            )
        )
    )

    fmt_kwargs = {"datefmt": datefmt}
    if style:
        vers = sys.version_info
        if vers < (3, 2):
            logging.warning(
                "Insufficient Python version to specify logging format style: "
                "{}.{}.{}".format(vers.major, vers.minor, vers.micro)
            )
        else:
            fmt_kwargs["style"] = style

    for h in handlers:
        h.setFormatter(logging.Formatter(get_fmt(h), **fmt_kwargs))
        h.setLevel(level)
        logger.addHandler(h)
    logger.debug(
        "Configured logger '%s' using %s v%s", logger.name, PACKAGE_NAME, __version__
    )

    return logger


def setup_logger(
    name="",
    level=None,
    stream=None,
    logfile=None,
    make_root=None,
    propagate=False,
    silent=False,
    devmode=False,
    verbosity=None,
    fmt=None,
    datefmt=None,
    plain_format=False,
    style=None,
):
    """Old alias for init_logger for backwards compatibility"""
    warnings.warn("Please use init_logger in place of setup_logger", DeprecationWarning)
    return init_logger(
        name,
        level,
        stream,
        logfile,
        make_root,
        propagate,
        silent,
        devmode,
        verbosity,
        fmt,
        datefmt,
        plain_format,
        style,
    )


def _level_from_verbosity(verbosity):
    """
    Translation of verbosity into logging level.

    Log message count monotonically increases in verbosity
    while it decreases in logging level, making verbosity
    a more intuitive specification mechanism for users.

    :param int | str verbosity: small integral value representing a relative
        measure of interest in seeing messages about program execution,
        or the name of a Python builtin logging level
    :return int: numeric logging level in accordance with Python builtin logging
    """
    try:
        verbosity = int(verbosity)
    except:
        pass
    if isinstance(verbosity, str):
        v = verbosity.upper()
        if v.startswith(_WARN_REPR):
            v = _WARN_REPR
        if v not in LEVEL_BY_VERBOSITY:
            raise ValueError(
                "Invalid logging verbosity ('{}'); choose from: "
                "{}".format(verbosity, ", ".join(LEVEL_BY_VERBOSITY))
            )
        return getattr(logging, v)
    elif isinstance(verbosity, int):
        return LEVEL_BY_VERBOSITY[verbosity - 1]  # 1-based user, 0-based internal
    else:
        raise TypeError(
            "Verbosity must be string or int; got {} ({})".format(
                verbosity, type(verbosity)
            )
        )


class AbsentOptionException(Exception):
    """Exception subtype suggesting that client should add log options."""

    def __init__(self, missing_optname):
        likely_reason = (
            "'{}' not in the parsed options; was {} used to "
            "add CLI logging options to an argument parser?".format(
                missing_optname, "{}.{}".format(__name__, add_logging_options.__name__)
            )
        )
        super(AbsentOptionException, self).__init__(likely_reason)


# Stolen from peppy. Probably need to make peppy/looper rely on this.
def get_logger(name):
    """
    Return a logger with given name, equipped with custom method.

    :param str name: name for the logger to get/create.
    :return logging.Logger: named, custom logger instance.
    """
    l = logging.getLogger(name)
    l.whisper = lambda msg, *args, **kwargs: l.log(5, msg, *args, **kwargs)
    return l
