import logging
from collections.abc import Mapping
from inspect import getfullargspec
from warnings import warn

from .const import *
from .exceptions import *
from .yacman import YacAttMap, _warn_deprecated

_LOGGER = logging.getLogger(__name__)


class AliasedYacAttMap(YacAttMap):
    """
    A class that extends YacAttMap to provide alias feature.

    The items in the object can be accessed using the original key or an alias,
    if defined in the aliases Mapping.
    """

    def __init__(
        self,
        entries=None,
        filepath=None,
        yamldata=None,
        writable=False,
        wait_max=DEFAULT_WAIT_TIME,
        skip_read_lock=False,
        aliases=None,
        exact=False,
        aliases_strict=None,
    ):
        """
        Object constructor

        :param Iterable[(str, object)] | Mapping[str, object] entries: YAML
            collection of key-value pairs.
        :param str filepath: YAML filepath to the config file.
        :param str yamldata: YAML-formatted string
        :param bool writable: whether to create the object with write
            capabilities
        :param int wait_max: how long to wait for creating an object when the
            file that data will be read from is locked
        :param bool skip_read_lock: whether the file should not be locked for
            reading when object is created in read only mode
        :param Mapping | callable(self) -> Mapping aliases: aliases mapping to
            use or a callable that produces such a mapping out of the object
            to set the aliases for
        :param bool aliases_strict: how to handle aliases mapping issues;
            None for warning, True for AliasError, False to disregard
        :param bool exact: whether aliases should not be used, even if defined
        """

        super(AliasedYacAttMap, self).__init__(
            entries=entries,
            filepath=filepath,
            yamldata=yamldata,
            writable=writable,
            wait_max=wait_max,
            skip_read_lock=skip_read_lock,
        )
        setattr(self[IK], ALIASES_KEY_RAW, {})
        if not exact:
            if isinstance(aliases, Mapping) and is_aliases_mapping_valid(
                aliases, aliases_strict
            ):
                setattr(self[IK], ALIASES_KEY_RAW, aliases)
            elif callable(aliases):
                if len(getfullargspec(aliases).args) != 1:
                    _emit_msg(
                        aliases_strict,
                        "Provided function '{}' must be a one-arg function".format(
                            aliases.__name__
                        ),
                    )
                try:
                    res = aliases(self)
                except Exception as e:
                    _emit_msg(
                        aliases_strict,
                        "Provided function '{}' errored: {}".format(
                            aliases.__name__, getattr(e, "message", repr(e))
                        ),
                    )
                else:
                    if is_aliases_mapping_valid(res):
                        setattr(self[IK], ALIASES_KEY_RAW, res)
                    else:
                        _emit_msg(
                            aliases_strict,
                            "callable '{}' did not return a Mapping".format(
                                aliases.__name__
                            ),
                        )
            else:
                _LOGGER.info("No aliases provided")

        # convert the original, condensed mapping to a data structure with
        # optimal time complexity
        setattr(self[IK], ALIASES_KEY, {})
        for k, v in self[IK][ALIASES_KEY_RAW].items():
            for alias in v:
                self[IK][ALIASES_KEY][alias] = k

    def __getitem__(self, item, expand=True, to_dict=False):
        """
        This item accession method will try to access the value by a literal
        key. If the key is not defined in the object it will try to access the
        key by it's alias, if defined. If both fail, a KeyError is raised.
        """
        try:
            return super(AliasedYacAttMap, self).__getitem__(
                item=item, expand=expand, to_dict=to_dict
            )
        except KeyError:
            try:
                key = self.get_key(item)
            except UndefinedAliasError:
                raise KeyError(item)
            else:
                return super(AliasedYacAttMap, self).__getitem__(
                    item=key, expand=expand, to_dict=to_dict
                )

    def __contains__(self, key):
        """
        This containment verification method will first try the  literal key.
        If the key is not defined in the object it will try to use its alias.
        If both fail, a negative decision is returned; otherwise -- positive.
        """
        try:
            self.__getitem__(key, expand=False)
        except (UndefinedAliasError, KeyError):
            try:
                alias_key = self.get_key(key)
            except (UndefinedAliasError, KeyError):
                return False
            else:
                try:
                    self.__getitem__(alias_key, expand=False)
                except (UndefinedAliasError, KeyError):
                    return False
                else:
                    return True
        else:
            return True

    def __delitem__(self, key):
        """
        This item deletion method will try to remove the item by the literal
        key or by its alias, if defined.
        """
        try:
            # check whether an alias was used
            alias_key = self.get_key(alias=key)
        except (UndefinedAliasError, KeyError):
            # alias was not used, try to delete the literal key
            super(AliasedYacAttMap, self).__delitem__(key=key)
        else:
            # alias was used, try to delete the alias
            super(AliasedYacAttMap, self).__delitem__(key=alias_key)

    @property
    def alias_dict(self):
        """
        Get the alias mapping bound to the object

        :return dict: alias-key mapping (one to one)
        """
        _warn_deprecated(obj=self)
        return getattr(self[IK], ALIASES_KEY)

    @property
    def _raw_alias_dict(self):
        """
        Get the raw alias mapping provided at the object construction.
        It is not updated, so it should not be used to the key-alias lookup

        :return dict: key-aliases mapping (one to many)
        """
        _warn_deprecated(obj=self)
        return getattr(self[IK], ALIASES_KEY_RAW)

    def get_aliases(self, key):
        """
        Get the alias for key in the object

        :param str key: key to find an alias for
        :return list[str]: aliases matched by the key
        :raise GenomeConfigFormatError: if aliases mapping has not been defined
            for this object
        :raise UndefinedAliasError: if no alias has been defined for the
            requested key
        """
        aliases = []
        for a, k in self[IK][ALIASES_KEY].items():
            if k == key:
                aliases.append(a)
        if aliases:
            return aliases
        raise UndefinedAliasError("No alias defined for: {}".format(key))

    def get_key(self, alias):
        """
        Get the key for alias in the object

        :param str alias: alias to find a key for
        :return str: key match by the alias
        :raise GenomeConfigFormatError: if aliases mapping has not been defined
            for this object
        :raise UndefinedAliasError: if a no key has been defined for the
            requested alias
        """
        try:
            # first try to use the parent method (doesn't try to use aliases) to
            # check if the __internal key is defined.
            # Otherwise we would end up in an infinite recursion loop.
            super(AliasedYacAttMap, self).__getitem__(item=IK, expand=False)
        except KeyError:
            # raise UndefinedAliasError, which is caught in the updated __getitem__ method
            raise UndefinedAliasError()
        if alias in self[IK][ALIASES_KEY].keys():
            return self[IK][ALIASES_KEY][alias]
        raise UndefinedAliasError("No key defined for: {}".format(alias))

    def set_aliases(self, key, aliases, overwrite=False, reset_key=False):
        """
        Assign an alias to a key in the object.

        :param str key: name of the key to assign to an alias for
        :param str | list[str] aliases: alias to use
        :param bool overwrite: whether to force overwrite the key for an
            already defined alias
        :param bool reset_key: whether to force remove existing aliases
            for a key
        :return list[str]: list of set aliases
        """
        removed_aliases = []
        if reset_key:
            try:
                current_aliases = self.get_aliases(key)
            except UndefinedAliasError:
                pass
            else:
                for a in current_aliases:
                    del self[IK][ALIASES_KEY][a]
                    removed_aliases.append(a)

        set_aliases = []
        for alias in _make_list_of_aliases(aliases):
            if alias in self[IK][ALIASES_KEY]:
                if overwrite:
                    self[IK][ALIASES_KEY][alias] = key
                    set_aliases.append(alias)
            else:
                self[IK][ALIASES_KEY][alias] = key
                set_aliases.append(alias)
        _LOGGER.debug("Added aliases ({}: {})".format(key, set_aliases))
        return set_aliases, removed_aliases

    def remove_aliases(self, key, aliases=None):
        """
        Remove an alias from the object.

        :param str key: name of the key to remove
        :param str aliases: list of aliases to remove
        :return list[str]: list of removed aliases
        """
        removed = []
        aliases = _make_list_of_aliases(aliases)
        try:
            current_aliases = self.get_aliases(key)
        except UndefinedAliasError:
            return removed
        else:
            existing_aliases = (
                list(set(aliases) & set(current_aliases))
                if aliases
                else current_aliases
            )
            for alias in existing_aliases:
                del self[IK][ALIASES_KEY][alias]
                removed.append(alias)
            return removed


def is_aliases_mapping_valid(aliases, strictness=None):
    """
    Determine if the aliases mapping is formatted properly, e.g. {"k": ["v"]}

    :param Mapping[list] aliases: mapping to verify
    :param bool strictness: how to handle format issues
    :return bool: whether the mapping adheres to the correct format
    """
    if isinstance(aliases, Mapping):
        if all([isinstance(v, list) for k, v in aliases.items()]):
            return True
    _emit_msg(strictness, "Provided aliases mapping is invalid; Mapping[list] required")
    return False


def _make_list_of_aliases(aliases):
    """
    Check and/or produce a proper aliases input

    :param str | list[str] aliases: alias or collection of aliases to check
    :return list[str]: list of aliases
    :raise AliasError: if the input format does not meet the requirements
    """

    def _raise_alias_class(x):
        raise AliasError(
            "A string or a list of strings is required, "
            "got: {}".format(x.__class__.__name__)
        )

    if aliases is None:
        return aliases
    if isinstance(aliases, str):
        aliases = [aliases]
    elif isinstance(aliases, list):
        assert all([isinstance(x, str) for x in aliases]), _raise_alias_class(aliases)
    else:
        _raise_alias_class(aliases)
    return aliases


def _emit_msg(strictness, msg):
    """
    Emit a message based on the selected strictness level

    :param bool strictness: None for warning, True for AliasError,
        False for debug log
    :param str msg: a message to emit
    """
    if strictness:
        raise AliasError(msg)
    elif strictness is None:
        warn(msg)
    else:
        _LOGGER.debug(msg)
