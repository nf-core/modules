""" AttMap that echoes an unset key/attr """

from .pathex_attmap import PathExAttMap

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"

__all__ = ["AttMapEcho", "EchoAttMap"]


class EchoAttMap(PathExAttMap):
    """An AttMap that returns key/attr if it has no set value."""

    def __getattr__(self, item, default=None, expand=True):
        """
        Fetch the value associated with the provided identifier.

        :param int | str item: identifier for value to fetch
        :param object default: default return value
        :param bool expand: whether to attempt variable expansion of string
            value, in case it's a path
        :return object: whatever value corresponds to the requested key/item
        :raises AttributeError: if the requested item has not been set,
            no default value is provided, and this instance is not configured
            to return the requested key/item itself when it's missing; also,
            if the requested item is unmapped and appears to be protected,
            i.e. by flanking double underscores, then raise AttributeError
            anyway. More specifically, respect attribute naming that appears
            to be indicative of the intent of protection.
        """
        try:
            return super(EchoAttMap, self).__getattr__(item, default, expand)
        except (AttributeError, TypeError):
            # If not, triage and cope accordingly.
            if self._is_od_member(item) or (
                item.startswith("__") and item.endswith("__")
            ):
                # Accommodate security-through-obscurity approach of some libs.
                error_reason = "Protected-looking attribute: {}".format(item)
                raise AttributeError(error_reason)
            return default if default is not None else item

    @property
    def _lower_type_bound(self):
        """Most specific type to which an inserted value may be converted"""
        return AttMapEcho


AttMapEcho = EchoAttMap
