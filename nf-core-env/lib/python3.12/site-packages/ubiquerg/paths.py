"""Filesystem utility functions"""

import os
import re

from typing import Any, Union, Optional

from .web import is_url

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


def expandpath(path):
    """Expand a filesystem path that may or may not contain user/env vars.

    Args:
        path: path to expand

    Returns:
        str: expanded version of input path
    """
    return os.path.expandvars(os.path.expanduser(path))


def parse_registry_path(
    rpstring: str,
    defaults: list[tuple[str, Any]] = [
        ("protocol", None),
        ("namespace", None),
        ("item", None),
        ("subitem", None),
        ("tag", None),
    ],
) -> Union[dict, None]:
    """Parse a 'registry path' string into components.

    A registry path is a string that is kind of like a URL, providing a unique
    identifier for a particular asset, like
    protocol::namespace/item.subitem:tag. You can use the `defaults` argument to
    change the names of the entries in the return dict, and to provide defaults
    in case of missing values.

    Args:
        rpstring: string to parse
        defaults: A list of 5 tuples with name of the 5 entries, and a default value in case it is missing (can be 'None')

    Returns:
        dict | None: dict with one element for each parsed entry in the path
    """

    # This commented regex is the same without protocol
    # ^(?:([0-9a-zA-Z_-]+)\/)?([0-9a-zA-Z_-]+)(?::([0-9a-zA-Z_.-]+))?$
    # regex = "^(?:([0-9a-zA-Z_-]+)(?:::|:\/\/))?(?:([0-9a-zA-Z_-]+)\/)?([0-9a-zA-Z_-]+)(?::([0-9a-zA-Z_.-]+))?$"
    regex = r"^(?:([0-9a-zA-Z._-]+)(?:::|:\/\/))?(?:([0-9a-zA-Z_-]+)\/)?([0-9a-zA-Z_-]+)(?:\.([0-9a-zA-Z_-]+))?(?::([0-9a-zA-Z_.,|+()-]+))?$"
    # This regex matches strings like:
    # protocol://namespace/item:tag
    # or: protocol::namespace/item:tag
    # The names 'protocol', 'namespace', 'item', and 'tag' are generic and
    # you can use this function for whatever you like in this format... The
    # regex can handle any of these missing and will parse correctly into the
    # same element
    # For instance, you can leave the tag or protocol or both off:
    # ucsc://hg38/bowtie2_index
    # hg38/bowtie2_index
    # With no delimiters, it will match the item name:
    # bowtie2_index

    res = re.match(regex, rpstring)
    if not res:
        return None
    # position 0: parent namespace
    # position 1: namespace
    # position 2: primary name
    # position 3: tag
    captures = res.groups()
    parsed_identifier = {
        defaults[0][0]: captures[0] or defaults[0][1],
        defaults[1][0]: captures[1] or defaults[1][1],
        defaults[2][0]: captures[2] or defaults[2][1],
        defaults[3][0]: captures[3] or defaults[3][1],
        defaults[4][0]: captures[4] or defaults[4][1],
    }
    return parsed_identifier


def parse_registry_path_strict(
    input_string: str,
    require_protocol: bool = False,
    require_namespace: bool = False,
    require_item: bool = True,
    require_subitem: bool = False,
    require_tag: bool = False,
) -> Union[dict[str, Any], None]:
    """Parse and validate a registry path with required component checks.

    This function parses a registry path and returns the parsed dictionary
    only if all required components are present. Returns None otherwise.
    Can be used as a boolean check (truthy/falsy) or to get the parsed components.

    Args:
        input_string: String to parse and validate as a registry path
        require_protocol: If True, protocol component must be present
        require_namespace: If True, namespace component must be present
        require_item: If True, item component must be present (default: True)
        require_subitem: If True, subitem component must be present
        require_tag: If True, tag component must be present

    Returns:
        dict | None: Parsed registry path dict if valid and all required components present, else None

    Example:
        >>> result = parse_registry_path_strict("namespace/item:tag")
        >>> result['namespace']
        'namespace'
        >>> parse_registry_path_strict("item", require_namespace=True)
        None
        >>> # Can be used as a boolean check
        >>> if parse_registry_path_strict("namespace/item", require_namespace=True):
        ...     print("Valid!")
        Valid!
        >>> # Get specific components
        >>> result = parse_registry_path_strict("protocol::namespace/item.subitem:tag", require_protocol=True)
        >>> result['protocol']
        'protocol'
    """
    parsed = parse_registry_path(input_string)

    if parsed is None:
        return None

    # Check required components
    requirements = {
        "protocol": require_protocol,
        "namespace": require_namespace,
        "item": require_item,
        "subitem": require_subitem,
        "tag": require_tag,
    }

    for component, required in requirements.items():
        if required and not parsed.get(component):
            return None

    return parsed


def mkabs(path: str, reldir: Optional[str] = None) -> str:
    """Make sure a path is absolute.

    If not already absolute, it's made absolute relative to a given directory (or file).
    Also expands ~ and environment variables for kicks.

    Args:
        path: Path to make absolute
        reldir: Relative directory to make path absolute from if it's not already absolute

    Returns:
        str: Absolute path
    """

    def xpand(path):
        return os.path.expandvars(os.path.expanduser(path))

    if path is None:
        return path

    if is_url(path):
        return path

    if os.path.isabs(xpand(path)):
        return xpand(path)

    if not reldir:
        return os.path.abspath(xpand(path))

    if os.path.isdir(reldir):
        return os.path.join(xpand(reldir), xpand(path))
    else:
        return os.path.join(xpand(os.path.dirname(reldir)), xpand(path))
