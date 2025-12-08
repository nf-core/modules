from logging import CRITICAL, DEBUG, ERROR, INFO, WARN

from peppy.const import SAMPLE_NAME_ATTR
from peppy import __version__ as peppy_version
from ubiquerg import VersionInHelpParser

from . import __version__
from .const import *

LEVEL_BY_VERBOSITY = [ERROR, CRITICAL, WARN, INFO, DEBUG]

version_combined = f"{__version__} (peppy {peppy_version})"


def build_argparser():
    banner = "%(prog)s - Interact with PEPs"
    additional_description = "\nhttp://eido.databio.org/"

    parser = VersionInHelpParser(
        prog=PKG_NAME,
        description=banner,
        epilog=additional_description,
        version=version_combined,
    )

    subparsers = parser.add_subparsers(dest="command")
    parser.add_argument(
        "--verbosity",
        dest="verbosity",
        type=int,
        choices=range(len(LEVEL_BY_VERBOSITY)),
        help="Choose level of verbosity (default: %(default)s)",
    )
    parser.add_argument("--logging-level", dest="logging_level", help="logging level")
    parser.add_argument(
        "--dbg",
        dest="dbg",
        action="store_true",
        help="Turn on debug mode (default: %(default)s)",
    )
    sps = {}
    for cmd, desc in SUBPARSER_MSGS.items():
        subparser = subparsers.add_parser(cmd, description=desc, help=desc)
        subparser.add_argument(
            "--st-index",
            required=False,
            type=str,
            default=SAMPLE_NAME_ATTR,
            help=f"Sample table index to use, samples are identified by '{SAMPLE_NAME_ATTR}' by default.",
        )
        subparser.add_argument(
            "--sst-index",
            required=False,
            type=str,
            default=SAMPLE_NAME_ATTR,
            help=f"Subsample table index to use, samples are identified by '{SAMPLE_NAME_ATTR}' by default.",
        )
        subparser.add_argument(
            "--amendments",
            required=False,
            type=str,
            nargs="+",
            help=f"Names of the amendments to activate.",
        )
        if cmd != CONVERT_CMD:
            subparser.add_argument(
                "pep",
                metavar="PEP",
                help="Path to a PEP configuration file in yaml format.",
                default=None,
            )
        else:
            subparser.add_argument(
                "pep",
                metavar="PEP",
                nargs="?",
                help="Path to a PEP configuration file in yaml format.",
                default=None,
            )

        sps[cmd] = subparser

    sps[VALIDATE_CMD].add_argument(
        "-s",
        "--schema",
        required=True,
        help="Path to a PEP schema file in yaml format.",
        metavar="S",
    )

    sps[INSPECT_CMD].add_argument(
        "-n",
        "--sample-name",
        required=False,
        nargs="+",
        help="Name of the samples to inspect.",
        metavar="SN",
    )

    sps[INSPECT_CMD].add_argument(
        "-l",
        "--attr-limit",
        required=False,
        type=int,
        default=10,
        help="Number of sample attributes to display.",
    )

    group = sps[VALIDATE_CMD].add_mutually_exclusive_group()

    group.add_argument(
        "-n",
        "--sample-name",
        required=False,
        help="Name or index of the sample to validate. "
        "Only this sample will be validated.",
        metavar="S",
    )

    group.add_argument(
        "-c",
        "--just-config",
        required=False,
        action="store_true",
        default=False,
        help="Whether samples should be excluded from the validation.",
    )

    sps[CONVERT_CMD].add_argument(
        "-f",
        "--format",
        required=False,
        default="yaml",
        help="Output format (name of filter; use -l to see available).",
    )

    sps[CONVERT_CMD].add_argument(
        "-n",
        "--sample-name",
        required=False,
        nargs="+",
        help="Name of the samples to inspect.",
    )

    sps[CONVERT_CMD].add_argument(
        "-a",
        "--args",
        nargs="+",
        action="append",
        required=False,
        default=None,
        help="Provide arguments to the filter function (e.g. arg1=val1 arg2=val2).",
    )

    sps[CONVERT_CMD].add_argument(
        "-l",
        "--list",
        required=False,
        default=False,
        action="store_true",
        help="List available filters.",
    )

    sps[CONVERT_CMD].add_argument(
        "-d",
        "--describe",
        required=False,
        default=False,
        action="store_true",
        help="Show description for a given filter.",
    )

    sps[CONVERT_CMD].add_argument(
        "-p",
        "--paths",
        nargs="+",
        help="Paths to dump conversion result as key=value pairs.",
    )
    return parser, sps
