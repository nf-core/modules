"""System utility functions"""

import os
import subprocess
from typing import Optional

__author__ = "Databio Lab"
__email__ = "nathan@code.databio.org"

__all__ = ["is_command_callable", "is_writable"]


def is_command_callable(cmd: str) -> bool:
    """Check if command can be called.

    Args:
        cmd: actual command to check for callability

    Returns:
        bool: whether given command's call succeeded

    Raises:
        TypeError: if the alleged command isn't a string
        ValueError: if the alleged command is empty
    """
    if not isinstance(cmd, str):
        raise TypeError(
            "Alleged command isn't a string: {} ({})".format(cmd, type(cmd))
        )
    if not cmd:
        raise ValueError("Empty command to check for callability")
    if os.path.isdir(cmd) or (os.path.isfile(cmd) and not os.access(cmd, os.X_OK)):
        return False
    # if system is windows run this command:
    if os.name == "nt":
        try:
            subprocess.run(
                ["where", cmd],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            return True
        except subprocess.CalledProcessError:
            return False
    else:
        # Use `command` to see if command is callable, and rule on exit code.
        check = "command -v {0} >/dev/null 2>&1 || {{ exit 1; }}".format(cmd)
        return not bool(os.system(check))


def is_writable(
    folder: Optional[str], check_exist: bool = False, create: bool = False
) -> bool:
    """Make sure a folder is writable.

    Given a folder, check that it exists and is writable. Errors if requested on
    a non-existent folder. Otherwise, make sure the first existing parent folder
    is writable such that this folder could be created.

    Args:
        folder: Folder to check for writeability
        check_exist: Throw an error if it doesn't exist?
        create: Create the folder if it doesn't exist?
    """
    folder = folder or "."

    if os.path.exists(folder):
        return os.access(folder, os.W_OK) and os.access(folder, os.X_OK)
    elif create:
        os.mkdir(folder)
        return True
    elif check_exist:
        raise OSError("Folder not found: {}".format(folder))
    else:
        # The folder didn't exist. Recurse up the folder hierarchy to make sure
        # all paths are writable
        return is_writable(os.path.dirname(folder), check_exist)
