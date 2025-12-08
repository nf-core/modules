"""Tools for working with collections"""

import itertools
import sys
from warnings import warn
from typing import Any, Optional, TypeVar
from collections.abc import Iterable, Mapping

T = TypeVar("T")

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


__all__ = [
    "is_collection_like",
    "powerset",
    "asciify_dict",
    "merge_dicts",
    "deep_update",
]


def merge_dicts(x: dict[Any, Any], y: dict[Any, Any]) -> dict[Any, Any]:
    """Merge dictionaries.

    Args:
        x: dict to merge
        y: dict to merge

    Returns:
        Mapping: merged dict
    """
    z = x.copy()
    z.update(y)
    return z


def deep_update(old: dict[Any, Any], new: Mapping[Any, Any]) -> dict[Any, Any]:
    """Recursively update nested dict, modifying source.

    Args:
        old: dict to update
        new: dict with new values

    Returns:
        dict: updated dict
    """
    for key, value in new.items():
        if isinstance(value, Mapping) and value:
            old[key] = deep_update(old.get(key, {}), value)
        else:
            old[key] = new[key]
    return old


def is_collection_like(c: Any) -> bool:
    """Determine whether an object is collection-like.

    Args:
        c: Object to test as collection

    Returns:
        bool: Whether the argument is a (non-string) collection
    """
    return isinstance(c, Iterable) and not isinstance(c, str)


def uniqify(seq: list[T]) -> list[T]:  # Dave Kirby
    """Return only unique items in a sequence, preserving order.

    Args:
        seq: List of items to uniqify

    Returns:
        list[object]: Original list with duplicates removed
    """
    # Order preserving
    seen = set()
    # Use list comprehension for speed
    return [x for x in seq if x not in seen and not seen.add(x)]  # type: ignore[func-returns-value]


def powerset(
    items: Iterable[T],
    min_items: Optional[int] = None,
    include_full_pop: bool = True,
    nonempty: bool = False,
) -> list[tuple[T, ...]]:
    """Build the powerset of a collection of items.

    Args:
        items: "Pool" of all items, the population for which to build the power set
        min_items: Minimum number of individuals from the population to allow in any given subset
        include_full_pop: Whether to include the full population in the powerset (default True to accord with genuine definition)
        nonempty: force each subset returned to be nonempty

    Returns:
        list[object]: Sequence of subsets of the population, in nondecreasing size order

    Raises:
        TypeError: if minimum item count is specified but is not an integer
        ValueError: if minimum item count is insufficient to guarantee nonempty subsets
    """
    if min_items is None:
        min_items = 1 if nonempty else 0
    else:
        if not isinstance(min_items, int):
            raise TypeError(
                "Min items count for each subset isn't an integer: "
                "{} ({})".format(min_items, type(min_items))
            )
        if nonempty and min_items < 1:
            raise ValueError(
                "When minimum item count is {}, nonempty subsets "
                "cannot be guaranteed.".format(min_items)
            )
    # Account for iterable burn possibility; besides, collection should be
    # relatively small if building the powerset.
    items = list(items)
    n = len(items)
    if n == 0 or n < min_items:
        return []
    max_items = len(items) + 1 if include_full_pop else len(items)
    return list(
        itertools.chain.from_iterable(
            itertools.combinations(items, k) for k in range(min_items, max_items)
        )
    )


def asciify_dict(data: dict[Any, Any]) -> dict[Any, Any]:
    """Legacy Python 2 function - now just returns input unchanged.

    This function was used to convert unicode strings to ASCII in Python 2.
    In Python 3+, strings are unicode by default, so this is a no-op.

    Reference: https://gist.github.com/chris-hailstorm/4989643

    Args:
        data: Dictionary to process

    Returns:
        dict: The same dictionary (unchanged in Python 3.9+)

    Note:
        Deprecated: This function is kept for backward compatibility but does nothing in Python 3.9+.
    """
    warn(
        "asciify_dict is deprecated and does nothing in Python 3.9+",
        DeprecationWarning,
        stacklevel=2,
    )
    return data
