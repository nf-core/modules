""" Canonical behavior for attmap in pepkit projects """

import sys

if sys.version_info < (3, 4):
    from collections import Mapping
else:
    from collections.abc import Mapping

from ubiquerg import expandpath

from .ordattmap import OrdAttMap

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


__all__ = ["PathExAttMap"]


class PathExAttMap(OrdAttMap):
    """Used in pepkit projects, with Mapping conversion and path expansion"""

    def __getattribute__(self, item, expand=True):
        res = super(PathExAttMap, self).__getattribute__(item)
        return _safely_expand(res) if expand else res

    def __getattr__(self, item, default=None, expand=True):
        """
        Get attribute, accessing stored key-value pairs as needed.

        :param str item: name of attribute/key
        :param object default: value to return if requested attr/key is missing
        :param bool expand: whether to attempt path expansion of string value
        :return object: value bound to requested name
        :raise AttributeError: if requested item is unavailable
        """
        try:
            v = super(PathExAttMap, self).__getattribute__(item)
        except AttributeError:
            try:
                return self.__getitem__(item, expand)
            except KeyError:
                # Requested item is unknown, but request was made via
                # __getitem__ syntax, not attribute-access syntax.
                raise AttributeError(item)
        else:
            return _safely_expand(v) if expand else v

    def __getitem__(self, item, expand=True, to_dict=False):
        """
        Fetch the value of given key.

        :param hashable item: key for which to fetch value
        :param bool expand: whether to expand string value as path
        :return object: value mapped to given key, if available
        :raise KeyError: if the requested key is unmapped.
        """
        v = super(PathExAttMap, self).__getitem__(item)
        return _safely_expand(v, to_dict) if expand else v

    def get(self, k, default=None, expand=True):
        try:
            return self.__getitem__(k, expand)
        except KeyError:
            return default

    def items(self, expand=False, to_dict=False):
        """
        Produce list of key-value pairs, optionally expanding paths.

        :param bool expand: whether to expand paths
        :return Iterable[object]: stored key-value pairs, optionally expanded
        """
        return [(k, self.__getitem__(k, expand, to_dict)) for k in self]

    def values(self, expand=False):
        """
        Produce list of values, optionally expanding paths.

        :param bool expand: whether to expand paths
        :return Iterable[object]: stored values, optionally expanded
        """
        return [self.__getitem__(k, expand) for k in self]

    def _data_for_repr(self, expand=False):
        """
        Hook for extracting the data used in the object's text representation.

        :param bool expand: whether to expand paths
        :return Iterable[(hashable, object)]: collection of key-value pairs
            to include in object's text representation
        """
        return filter(
            lambda kv: not self._excl_from_repr(kv[0], self.__class__),
            self.items(expand),
        )

    def to_map(self, expand=False):
        """
        Convert this instance to a dict.

        :return dict[str, object]: this map's data, in a simpler container
        """
        return self._simplify_keyvalue(self.items(expand), self._new_empty_basic_map)

    def to_dict(self, expand=False):
        """
        Return a builtin dict representation of this instance.

        :return dict: builtin dict representation of this instance
        """
        return self._simplify_keyvalue(self.items(expand, to_dict=True), dict)

    @property
    def _lower_type_bound(self):
        return PathExAttMap


def _safely_expand(x, to_dict=False):
    if isinstance(x, str):
        return expandpath(x)
    if to_dict and isinstance(x, Mapping):
        return {k: _safely_expand(v, to_dict) for k, v in x.items()}
    return x
