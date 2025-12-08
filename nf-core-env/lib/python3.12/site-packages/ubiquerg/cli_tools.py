"""Functions for working with command-line interaction"""

from .collection import is_collection_like, merge_dicts
from argparse import (
    ArgumentParser,
    _SubParsersAction,
    _HelpAction,
    _VersionAction,
    SUPPRESS,
)
import sys
from typing import Any, Optional, Union, Iterable

__classes__ = ["VersionInHelpParser"]
__all__ = __classes__ + ["build_cli_extra", "query_yes_no", "convert_value"]


class VersionInHelpParser(ArgumentParser):
    def __init__(self, version: Optional[str] = None, **kwargs: Any) -> None:
        """Overwrites the inherited init. Saves the version as an object attribute for further use."""
        super(VersionInHelpParser, self).__init__(**kwargs)
        self.version = version
        if self.version is not None:
            self.add_argument(
                "--version",
                action="version",
                version="%(prog)s {}".format(self.version),
            )

    def format_help(self) -> str:
        """Add version information to help text."""
        help_string = (
            "version: {}\n".format(str(self.version))
            if self.version is not None
            else ""
        )
        return help_string + super(VersionInHelpParser, self).format_help()

    def subparsers(self) -> _SubParsersAction:
        """Get the subparser associated with a parser.

        Returns:
            argparse._SubparsersAction: action defining the subparsers
        """
        subs = [a for a in self._actions if isinstance(a, _SubParsersAction)]
        if len(subs) != 1:
            raise ValueError("Expected exactly 1 subparser, got {}".format(len(subs)))
        return subs[0]

    def top_level_args(self) -> list[Any]:
        """Get actions not assiated with any subparser.

        Help and version are also excluded.

        Returns:
            list[argparse.<action_type>]: list of argument actions
        """
        excl = [_SubParsersAction, _HelpAction, _VersionAction]
        return [a for a in self._actions if not type(a) in excl]

    def subcommands(self) -> list[str]:
        """Get subcommands defined by a parser.

        Returns:
            list[str]: subcommands defined within this parser
        """
        return list(self.subparsers().choices.keys())

    def dests_by_subparser(
        self, subcommand: Optional[str] = None, top_level: bool = False
    ) -> Union[list[str], dict[str, list[str]]]:
        """Get argument dests by subcommand from a parser.

        Args:
            subcommand: subcommand to get dests for

        Returns:
            dict: dests by subcommand
        """
        if top_level:
            top_level_actions = self.top_level_args()
            dest_list = []
            for tla in top_level_actions:
                if hasattr(tla, "dest"):
                    dest_list.append(tla.dest)
            return dest_list

        if subcommand is not None and subcommand not in self.subcommands():
            raise ValueError(
                "'{}' not found in this parser commands: {}".format(
                    subcommand, str(self.subcommands())
                )
            )
        subs = (
            self.subparsers().choices
            if subcommand is None
            else {subcommand: self.subparsers().choices[subcommand]}
        )
        dests = {}
        for subcmd, sub in subs.items():
            dest_list: list[str] = []
            for action in sub._actions:
                if isinstance(action, _HelpAction):
                    continue
                if hasattr(action, "dest"):
                    dest_list.append(action.dest)
            dests[subcmd] = dest_list
        return dests

    def suppress_defaults(self) -> None:
        """Remove parser change defaults to argparse.SUPPRESS.

        This prevents them from showing up in the argparse.Namespace object after argument parsing.
        """
        top_level_actions = self.top_level_args()
        for tla in top_level_actions:
            if hasattr(tla, "dest"):
                tla.dest = SUPPRESS
        subs = self.subparsers().choices
        for subcmd, sub in subs.items():
            for sa in sub._actions:
                if hasattr(sa, "dest"):
                    sa.default = SUPPRESS

    def arg_defaults(
        self,
        subcommand: Optional[str] = None,
        unique: bool = False,
        top_level: bool = False,
    ) -> Union[dict[str, Any], dict[str, dict[str, Any]]]:
        """Get argument defaults by subcommand from a parser.

        Args:
            subcommand: subcommand to get defaults for
            unique: whether only unique flat dict of dests and defaults mapping should be returned

        Returns:
            dict: defaults by subcommand
        """
        if top_level:
            top_level_actions = self.top_level_args()
            defaults_dict = {}
            for tla in top_level_actions:
                if hasattr(tla, "default") and hasattr(tla, "dest"):
                    defaults_dict.update({tla.dest: tla.default})
            return defaults_dict

        if subcommand is not None and subcommand not in self.subcommands():
            raise ValueError(
                "'{}' not found in this parser commands: {}".format(
                    subcommand, str(self.subcommands())
                )
            )
        subs = (
            self.subparsers().choices
            if subcommand is None
            else {subcommand: self.subparsers().choices[subcommand]}
        )
        defaults = {}
        for subcmd, sub in subs.items():
            defaults_dict = {}
            for action in sub._actions:
                if isinstance(action, _HelpAction):
                    continue
                if hasattr(action, "default") and hasattr(action, "dest"):
                    defaults_dict.update({action.dest: action.default})
            defaults[subcmd] = defaults_dict
            if unique:
                unique_defaults = {}
                for k, v in defaults.items():
                    unique_defaults = merge_dicts(unique_defaults, v)
                return unique_defaults
        return defaults


def build_cli_extra(optargs: Union[dict[str, Any], Iterable]) -> str:
    """Render CLI options/args as text to add to base command.

    To specify a flag, map an option to None. Otherwise, map option short or
    long name to value(s). Values that are collection types will be rendered
    with single space between each. All non-string values are converted to
    string.

    Args:
        optargs: values used as options/arguments

    Returns:
        str: text to add to base command, based on given opts/args

    Raises:
        TypeError: if an option name isn't a string
    """

    def render(k, v):
        if not isinstance(k, str):
            raise TypeError("Option name isn't a string: {} ({})".format(k, type(k)))
        if v is None:
            return k
        if is_collection_like(v):
            v = " ".join(map(str, v))
        return "{} {}".format(k, v)

    try:
        data_iter = optargs.items()
    except AttributeError:
        data_iter = optargs

    return " ".join(render(*kv) for kv in data_iter)


def query_yes_no(question: str, default: str = "no") -> bool:
    """Ask a yes/no question via input() and return their answer.

    Args:
        question: a string that is presented to the user.
        default: the presumed answer if the user just hits <Enter>.

    Returns:
        bool: True for "yes" or False for "no"
    """

    def parse(ans):
        return {"yes": True, "y": True, "ye": True, "no": False, "n": False}[
            ans.lower()
        ]

    try:
        prompt = {None: "[y/n]", "yes": "[Y/n]", "no": "[y/N]"}[
            None if default is None else default.lower()
        ]
    except (AttributeError, KeyError):
        raise ValueError("invalid default answer: {}".format(default))
    msg = "{q} {p} ".format(q=question, p=prompt)
    while True:
        sys.stdout.write(msg)
        try:
            return parse(input() or default)
        except KeyError:
            sys.stdout.write("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")


def convert_value(val: Any) -> Union[bool, str, int, float, None]:
    """Convert string to the most appropriate type.

    Converts to one of: bool, str, int, None or float

    Args:
        val: the string to convert

    Returns:
        bool | str | int | float | None: converted string to the most appropriate type
    """
    if not isinstance(val, str):
        try:
            val = str(val)
        except:
            raise ValueError(
                "The input has to be of type convertible to 'str',"
                " got '{}'".format(type(val))
            )

    # val is definitely a string at this point
    if val == "None":
        return None
    if val.lower() == "true":
        return True
    if val.lower() == "false":
        return False

    try:
        return int(val)
    except ValueError:
        try:
            return float(val)
        except ValueError:
            return val
