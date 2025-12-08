"""Web-related utilities"""

import re

__all__ = ["is_url"]


def is_url(maybe_url: str) -> bool:
    """Determine whether a path is a URL.

    Args:
        maybe_url: path to investigate as URL

    Returns:
        bool: whether path appears to be a URL
    """
    # from Django 1.3.x
    # https://github.com/django/django/blob/6726d750979a7c29e0dd866b4ea367eef7c8a420/django/core/validators.py#L45-L51
    regex = re.compile(
        r"^(?:http|ftp)s?://"
        r"(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|"
        r"localhost|"
        r"\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})"
        r"(?::\d+)?"
        r"(?:/?|[/?]\S+)$",
        re.IGNORECASE,
    )
    return re.match(regex, str(maybe_url)) is not None
