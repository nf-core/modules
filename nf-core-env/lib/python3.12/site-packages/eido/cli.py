import logging
import sys

from logmuse import init_logger
from peppy import Project

from .argparser import LEVEL_BY_VERBOSITY, build_argparser
from .const import *
from .conversion import (
    convert_project,
    get_available_pep_filters,
    pep_conversion_plugins,
)
from .exceptions import EidoFilterError, EidoValidationError
from .inspection import inspect_project
from .validation import validate_config, validate_project, validate_sample


def _parse_filter_args_str(input):
    """
    Parse user input specification.

    :param Iterable[Iterable[str]] input: user command line input,
        formatted as follows: [[arg=txt, arg1=txt]]
    :return dict: mapping of keys, which are input names and values
    """
    lst = []
    for i in input or []:
        lst.extend(i)
    return (
        {x.split("=")[0]: x.split("=")[1] for x in lst if "=" in x}
        if lst is not None
        else lst
    )


def print_error_summary(errors_by_type):
    """Print a summary of errors, organized by error type"""
    n_error_types = len(errors_by_type)
    print(f"Found {n_error_types} types of error:")
    for type in errors_by_type:
        n = len(errors_by_type[type])
        msg = f"  - {type}: ({n} samples) "
        if n < 50:
            msg += ", ".join([x["sample_name"] for x in errors_by_type[type]])
        print(msg)

    if len(errors_by_type) > 1:
        final_msg = f"Validation unsuccessful. {len(errors_by_type)} error types found."
    else:
        final_msg = f"Validation unsuccessful. {len(errors_by_type)} error type found."

    print(final_msg)
    return final_msg


def main():
    """Primary workflow"""
    parser, sps = build_argparser()
    args, remaining_args = parser.parse_known_args()

    if args.command is None:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Set the logging level.
    if args.dbg:
        # Debug mode takes precedence and will listen for all messages.
        level = args.logging_level or logging.DEBUG
    elif args.verbosity is not None:
        # Verbosity-framed specification trumps logging_level.
        level = LEVEL_BY_VERBOSITY[args.verbosity]
    else:
        # Normally, we're not in debug mode, and there's not verbosity.
        level = LOGGING_LEVEL

    logger_kwargs = {"level": level, "devmode": args.dbg}
    init_logger(name="peppy", **logger_kwargs)
    global _LOGGER
    _LOGGER = init_logger(name=PKG_NAME, **logger_kwargs)

    if args.command == CONVERT_CMD:
        filters = get_available_pep_filters()
        if args.list:
            _LOGGER.info("Available filters:")
            if len(filters) < 1:
                _LOGGER.info("No available filters")
            for filter_name in filters:
                _LOGGER.info(f" - {filter_name}")
            sys.exit(0)
        if not "format" in args:
            _LOGGER.info("The following arguments are required: --format")
            sps[CONVERT_CMD].print_help(sys.stderr)
            sys.exit(1)
        if args.describe:
            if args.format not in filters:
                raise EidoFilterError(
                    f"'{args.format}' filter not found. Available filters: {', '.join(filters)}"
                )
            filter_functions_by_name = pep_conversion_plugins()
            print(filter_functions_by_name[args.format].__doc__)
            sys.exit(0)
        if args.pep is None:
            sps[CONVERT_CMD].print_help(sys.stderr)
            _LOGGER.info("The following arguments are required: PEP")
            sys.exit(1)
        if args.paths:
            paths = {y[0]: y[1] for y in [x.split("=") for x in args.paths]}
        else:
            paths = None

        p = Project(
            args.pep,
            sample_table_index=args.st_index,
            subsample_table_index=args.sst_index,
            amendments=args.amendments,
        )
        plugin_kwargs = _parse_filter_args_str(args.args)

        # append paths
        plugin_kwargs["paths"] = paths

        convert_project(p, args.format, plugin_kwargs)
        _LOGGER.info("Conversion successful")
        sys.exit(0)

    _LOGGER.debug(f"Creating a Project object from: {args.pep}")
    if args.command == VALIDATE_CMD:
        p = Project(
            args.pep,
            sample_table_index=args.st_index,
            subsample_table_index=args.sst_index,
            amendments=args.amendments,
        )
        if args.sample_name:
            try:
                args.sample_name = int(args.sample_name)
            except ValueError:
                pass
            _LOGGER.debug(
                f"Comparing Sample ('{args.pep}') in Project ('{args.pep}') "
                f"against a schema: {args.schema}"
            )
            validator = validate_sample
            arguments = [p, args.sample_name, args.schema]
        elif args.just_config:
            _LOGGER.debug(
                f"Comparing Project ('{args.pep}') against a schema: {args.schema}"
            )
            validator = validate_config
            arguments = [p, args.schema]
        else:
            _LOGGER.debug(
                f"Comparing Project ('{args.pep}') against a schema: {args.schema}"
            )
            validator = validate_project
            arguments = [p, args.schema]
        try:
            validator(*arguments)
        except EidoValidationError as e:
            print_error_summary(e.errors_by_type)
            return False
        _LOGGER.info("Validation successful")
        sys.exit(0)

    if args.command == INSPECT_CMD:
        p = Project(
            args.pep,
            sample_table_index=args.st_index,
            subsample_table_index=args.sst_index,
            amendments=args.amendments,
        )
        inspect_project(p, args.sample_name, args.attr_limit)
        sys.exit(0)
