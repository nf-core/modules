""" The trait defining a multi-access data object """

import abc
import sys

if sys.version_info < (3, 3):
    from collections import Mapping, MutableMapping
else:
    from collections.abc import Mapping, MutableMapping

from .helpers import get_data_lines, get_logger, is_custom_map

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"

__all__ = ["AttMapLike"]


_LOGGER = get_logger(__name__)


class AttMapLike(MutableMapping):
    """Base class for multi-access-mode data objects."""

    __metaclass__ = abc.ABCMeta

    def __init__(self, entries=None):
        """
        Create a new instance, optionally with initial key-value pairs.

        :param Mapping | Iterable[(Hashable, object)] entries: initial
            KV pairs to store
        """
        self.add_entries(entries)

    def __getattr__(self, item, default=None):
        try:
            return super(AttMapLike, self).__getattribute__(item)
        except AttributeError:
            try:
                return self.__getitem__(item)
            except KeyError:
                # Requested item is unknown, but request was made via
                # __getitem__ syntax, not attribute-access syntax.
                raise AttributeError(item)

    @abc.abstractmethod
    def __delitem__(self, item):
        pass

    @abc.abstractmethod
    def __getitem__(self, item):
        pass

    @abc.abstractmethod
    def __setitem__(self, key, value):
        pass

    def __iter__(self):
        return iter([k for k in self.__dict__.keys()])

    def __len__(self):
        return sum(1 for _ in iter(self))

    def __repr__(self):
        return self._render(
            self._simplify_keyvalue(self._data_for_repr(), self._new_empty_basic_map)
        )

    def _render(self, data, exclude_class_list=[]):
        def _custom_repr(obj, prefix=""):
            """
            Calls the ordinary repr on every object but list, which is
            converted to a block style string instead.

            :param object obj: object to convert to string representation
            :param str prefix: string to prepend to each list line in block
            :return str: custom object representation
            """
            if isinstance(obj, list) and len(obj) > 0:
                return f"\n{prefix} - " + f"\n{prefix} - ".join([str(i) for i in obj])
            return obj.strip("'") if hasattr(obj, "strip") else str(obj)

        class_name = self.__class__.__name__
        if class_name in exclude_class_list:
            base = ""
        else:
            base = class_name + "\n"

        if data:
            return base + "\n".join(get_data_lines(data, _custom_repr))
        else:
            return class_name + ": {}"

    def add_entries(self, entries):
        """
        Update this instance with provided key-value pairs.

        :param Iterable[(object, object)] | Mapping | pandas.Series entries:
            collection of pairs of keys and values
        """
        if entries is None:
            return
        # Permit mapping-likes and iterables/generators of pairs.
        if callable(entries):
            entries = entries()
        elif any("pandas.core" in str(t) for t in type(entries).__bases__):
            entries = entries.to_dict()
        try:
            entries_iter = entries.items()
        except AttributeError:
            entries_iter = entries
        for k, v in entries_iter:
            self[k] = (
                v
                if (
                    k not in self
                    or not isinstance(v, Mapping)
                    or not isinstance(self[k], Mapping)
                )
                else self[k].add_entries(v)
            )
        return self

    def get_yaml_lines(
        self,
        conversions=((lambda obj: isinstance(obj, Mapping) and 0 == len(obj), None),),
    ):
        """
        Get collection of lines that define YAML text rep. of this instance.

        :param Iterable[(function(object) -> bool, object)] conversions:
            collection of pairs in which first component is predicate function
            and second is what to replace a value with if it satisfies the predicate
        :return list[str]: YAML representation lines
        """
        if 0 == len(self):
            return ["{}"]
        data = self._simplify_keyvalue(
            self._data_for_repr(), self._new_empty_basic_map, conversions=conversions
        )
        return self._render(data).split("\n")[1:]

    def is_null(self, item):
        """
        Conjunction of presence in underlying mapping and value being None

        :param object item: Key to check for presence and null value
        :return bool: True iff the item is present and has null value
        """
        return item in self and self[item] is None

    def non_null(self, item):
        """
        Conjunction of presence in underlying mapping and value not being None

        :param object item: Key to check for presence and non-null value
        :return bool: True iff the item is present and has non-null value
        """
        return item in self and self[item] is not None

    def to_map(self):
        """
        Convert this instance to a dict.

        :return dict[str, object]: this map's data, in a simpler container
        """
        return self._simplify_keyvalue(self.items(), self._new_empty_basic_map)

    def to_dict(self):
        """
        Return a builtin dict representation of this instance.

        :return dict: builtin dict representation of this instance
        """
        return self._simplify_keyvalue(self.items(), dict)

    def to_yaml(self, trailing_newline=True):
        """
        Get text for YAML representation.

        :param bool trailing_newline: whether to add trailing newline
        :return str: YAML text representation of this instance.
        """
        return "\n".join(self.get_yaml_lines()) + ("\n" if trailing_newline else "")

    def _data_for_repr(self):
        """
        Hook for extracting the data used in the object's text representation.

        :return Iterable[(hashable, object)]: collection of key-value pairs
            to include in object's text representation
        """
        return filter(
            lambda kv: not self._excl_from_repr(kv[0], self.__class__), self.items()
        )

    def _excl_from_eq(self, k):
        """
        Hook for exclusion of particular value from a representation

        :param hashable k: key to consider for omission
        :return bool: whether the given key k should be omitted from comparison
        """
        return False

    def _excl_from_repr(self, k, cls):
        """
        Hook for exclusion of particular value from a representation

        :param hashable k: key to consider for omission
        :param type cls: data type on which to base the exclusion
        :return bool: whether the given key k should be omitted from
            text representation
        """
        return False

    def _excl_classes_from_todict(self):
        """
        Hook for exclusion of particular class from a dict conversion
        """
        return

    @abc.abstractproperty
    def _lower_type_bound(self):
        """Most specific type to which stored Mapping should be transformed"""
        pass

    @abc.abstractmethod
    def _new_empty_basic_map(self):
        """Return the empty collection builder for Mapping type simplification."""
        pass

    def _simplify_keyvalue(
        self,
        kvs,
        build,
        acc=None,
        conversions=None,
    ):
        """
        Simplify a collection of key-value pairs, "reducing" to simpler types.

        :param Iterable[(object, object)] kvs: collection of key-value pairs
        :param callable build: how to build an empty collection
        :param Iterable acc: accumulating collection of simplified data
        :return Iterable: collection of simplified data
        """
        acc = acc or build()
        kvs = iter(kvs)
        try:
            k, v = next(kvs)
        except StopIteration:
            return acc
        if not isinstance(v, self._excl_classes_from_todict() or tuple()):
            if is_custom_map(v):
                v = self._simplify_keyvalue(v.items(), build, build())
            if isinstance(v, Mapping):
                for pred, proxy in conversions or []:
                    if pred(v):
                        v = proxy
                        break
            acc[k] = v
        return self._simplify_keyvalue(kvs, build, acc, conversions)
