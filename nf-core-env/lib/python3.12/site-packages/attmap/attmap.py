""" Dot notation support for Mappings. """

import sys

if sys.version_info < (3, 3):
    from collections import Mapping
else:
    from collections.abc import Mapping

from ._att_map_like import AttMapLike
from .helpers import copy, get_logger, safedel_message

_LOGGER = get_logger(__name__)


@copy
class AttMap(AttMapLike):
    """
    A class to convert a nested mapping(s) into an object(s) with key-values
    using object syntax (attmap.attribute) instead of getitem syntax
    (attmap["key"]). This class recursively sets mappings to objects,
    facilitating attribute traversal (e.g., attmap.attr.attr).
    """

    def __delitem__(self, key):
        try:
            del self.__dict__[key]
        except KeyError:
            _LOGGER.debug(safedel_message(key))

    def __getitem__(self, item):
        return self.__dict__[item]

    def __setitem__(self, key, value):
        """
        This is the key to making this a unique data type.

        :param str key: name of the key/attribute for which to establish value
        :param object value: value to which set the given key; if the value is
            a mapping-like object, other keys' values may be combined.
        """
        # TODO: consider enforcement of type constraint, that value of different
        # type may not overwrite existing.
        self.__dict__[key] = self._final_for_store(key, value)

    def __eq__(self, other):
        # TODO: check for equality across classes?
        if (type(self) != type(other)) or (len(self) != len(other)):
            return False
        for k, v in self.items():
            if self._excl_from_eq(k):
                _LOGGER.debug("Excluding from comparison: {}".format(k))
                continue
            if not self._cmp(v, other[k]):
                return False
        return True

    def __ne__(self, other):
        return not self == other

    @staticmethod
    def _cmp(a, b):
        """Hook to tailor value comparison in determination of map equality."""

        def same_type(obj1, obj2, typenames=None):
            t1, t2 = str(obj1.__class__), str(obj2.__class__)
            return (t1 in typenames and t2 in typenames) if typenames else t1 == t2

        if same_type(
            a, b, ["<type 'numpy.ndarray'>", "<class 'numpy.ndarray'>"]
        ) or same_type(a, b, ["<class 'pandas.core.series.Series'>"]):
            check = lambda x, y: (x == y).all()
        elif same_type(a, b, ["<class 'pandas.core.frame.DataFrame'>"]):
            check = lambda x, y: (x == y).all().all()
        else:
            check = lambda x, y: x == y
        try:
            return check(a, b)
        except ValueError:
            # ValueError arises if, e.g., the pair of Series have
            # have nonidentical labels.
            return False

    def _final_for_store(self, k, v):
        """
        Before storing a value, apply any desired transformation.

        :param hashable k: key for which to store value
        :param object v: value to potentially transform before storing
        :return object: finalized value
        """
        if isinstance(v, Mapping) and not isinstance(v, self._lower_type_bound):
            v = self._metamorph_maplike(v)
        return v

    @property
    def _lower_type_bound(self):
        return AttMap

    def _metamorph_maplike(self, m):
        """
        Ensure a stored Mapping conforms with type expectation.

        :param Mapping m: the mapping to which to apply type transformation
        :return Mapping: a (perhaps more specialized) version of the given map
        :raise TypeError: if the given value isn't a Mapping
        """
        if not isinstance(m, Mapping):
            raise TypeError(
                "Cannot integrate a non-Mapping: {}\nType: {}".format(m, type(m))
            )
        return self._lower_type_bound(m.items())

    def _new_empty_basic_map(self):
        """Return the empty collection builder for Mapping type simplification."""
        return dict()

    def _repr_pretty_(self, p, cycle):
        """
        IPython display; https://ipython.readthedocs.io/en/stable/api/generated/IPython.lib.pretty.html

        :param IPython.lib.pretty.PrettyPrinter p: printer instance
        :param bool cycle: whether a cyclic reference is detected
        :return str: text representation of the instance
        """
        return p.text(repr(self) if not cycle else "...")
