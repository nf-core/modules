import logging
import os
import warnings
from collections.abc import Iterable, Mapping
from sys import _getframe
from urllib.request import urlopen
from warnings import warn

import attmap
import oyaml as yaml
from jsonschema import validate as _validate
from jsonschema.exceptions import ValidationError
from ubiquerg import create_lock, expandpath, is_url, make_lock_path, mkabs, remove_lock

from .const import *
from ._version import __version__
from typing import Union
from pathlib import Path


_LOGGER = logging.getLogger(__name__)

# Hack for string indexes of both ordered and unordered yaml representations
# Credit: Anthon
# https://stackoverflow.com/questions/50045617
# https://stackoverflow.com/questions/5121931
# The idea is: if you have yaml keys that can be interpreted as an int or a float,
# then the yaml loader will convert them into an int or a float, and you would
# need to access them with dict[2] instead of dict['2']. But since we always
# expect the keys to be strings, this doesn't work. So, here we are adjusting
# the loader to keep everything as a string. This happens in 2 ways, so that
# it's compatible with both yaml and oyaml, which is the orderedDict version.
# this will go away in python 3.7, because the dict representations will be
# ordered by default.

# Only do once?
if not hasattr(yaml.SafeLoader, "patched_yaml_loader"):
    _LOGGER.debug("Patching yaml loader")

    def my_construct_mapping(self, node, deep=False):
        data = self.construct_mapping_org(node, deep)
        return {
            (str(key) if isinstance(key, float) or isinstance(key, int) else key): data[
                key
            ]
            for key in data
        }

    yaml.SafeLoader.construct_mapping_org = yaml.SafeLoader.construct_mapping
    yaml.SafeLoader.construct_mapping = my_construct_mapping
    yaml.SafeLoader.patched_yaml_loader = True

# End hack


# from typing_extensions import deprecated
# @deprecated("YacAttMap is deprecated. Use YAMLConfigManager instead.")
class YacAttMap(attmap.PathExAttMap):
    """
    A class that extends AttMap to provide yaml reading and race-free
    writing in multi-user contexts.

    The YacAttMap class is a YAML Configuration Attribute Map. Think of it as a
    Python representation of your YAML configuration file, that can do a lot of cool
    stuff. You can access the hierarchical YAML attributes with dot notation or dict
    notation. You can read and write YAML config files with easy functions. It also
    retains memory of the its source filepath. If both a filepath and an entries
    dict are provided, it will first load the file and then updated it with
    values from the dict. Moreover, the config contents can be validated against a
    jsonschema schema, if a path to one is provided.
    """

    def __init__(
        self,
        entries=None,
        filepath=None,
        yamldata=None,
        writable=False,
        wait_max=DEFAULT_WAIT_TIME,
        skip_read_lock=False,
        schema_source=None,
        write_validate=False,
    ):
        """
        Object constructor

        :param Iterable[(str, object)] | Mapping[str, object] entries: YAML collection
            of key-value pairs.
        :param str filepath: YAML filepath to the config file.
        :param str yamldata: YAML-formatted string
        :param bool writable: whether to create the object with write capabilities
        :param int wait_max: how long to wait for creating an object when the file
            that data will be read from is locked
        :param bool skip_read_lock: whether the file should not be locked for reading
            when object is created in read only mode
        :param str schema_source: path or a URL to a jsonschema in YAML format to use
            for optional config validation. If this argument is provided the object
            is always validated at least once, at the object creation stage.
        :param bool write_validate: a boolean indicating whether the object should be
            validated every time the `write` method is executed, which is
            a way of preventing invalid config writing
        """
        if writable:
            if filepath:
                create_lock(filepath, wait_max)
            else:
                warnings.warn(
                    "Argument 'writable' is disregarded when the object is created "
                    "with 'entries' rather than read from the 'filepath'",
                    UserWarning,
                )
        if filepath:
            if not skip_read_lock and not writable and os.path.exists(filepath):
                create_lock(filepath, wait_max)
                file_contents = load_yaml(filepath)
                remove_lock(filepath)
            else:
                file_contents = load_yaml(filepath)
            if entries:
                if file_contents is None:
                    # if file is empty, initialize its contents to an empty dict
                    file_contents = {}
                file_contents.update(entries)
            entries = file_contents
        elif yamldata:
            entries = yaml.load(yamldata, yaml.SafeLoader)
        if not hasattr(self, IK):
            setattr(self, IK, attmap.AttMap())
        super(YacAttMap, self).__init__(entries or {})
        if filepath:
            # to make this python2 compatible, the attributes need to be set here.
            # prevents: AttributeError: _OrderedDict__root
            setattr(self[IK], WAIT_MAX_KEY, wait_max)
            setattr(self[IK], FILEPATH_KEY, mkabs(filepath))
            setattr(self[IK], RO_KEY, not writable)

        setattr(self[IK], WRITE_VALIDATE_KEY, write_validate)
        if schema_source is not None:
            assert isinstance(schema_source, str), TypeError(
                f"Path to the schema to validate the config must be a string"
            )
            sp = expandpath(schema_source)
            assert os.path.exists(sp), FileNotFoundError(
                f"Provided schema file does not exist: {schema_source}."
                f" Also tried: {sp}"
            )
            # validate config
            setattr(self[IK], SCHEMA_KEY, load_yaml(sp))
            self.validate()

    def __del__(self):
        try:
            if hasattr(self[IK], FILEPATH_KEY) and not getattr(self[IK], RO_KEY, True):
                self.make_readonly()
        except KeyError:
            # The __internal key may not exist during garbage collection
            pass

    def __repr__(self):
        # Here we want to render the data in a nice way; and we want to indicate
        # the class if it's NOT a YacAttMap. If it is a YacAttMap we just want
        # to give you the data without the class name.
        return self._render(
            self._simplify_keyvalue(self._data_for_repr(), self._new_empty_basic_map),
            exclude_class_list="YacAttMap",
        )

    def __enter__(self):
        setattr(self[IK], ORI_STATE_KEY, getattr(self[IK], RO_KEY, True))
        if not getattr(self[IK], RO_KEY, None):
            return self
        else:
            self.make_writable()
            return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.write()
        if getattr(self[IK], ORI_STATE_KEY, False):
            self.make_readonly()

    def _reinit(self, filepath=None):
        """
        Re-initialize the object

        :param str filepath: path to the file that should be read
        """
        if filepath is not None:
            self.__init__(filepath=filepath, skip_read_lock=True)
        else:
            self.__init__(entries={}, skip_read_lock=True)

    def _excl_from_repr(self, k, cls):
        return k in ATTR_KEYS

    @property
    def _lower_type_bound(self):
        """Most specific type to which an inserted value may be converted"""
        return YacAttMap

    def validate(self, schema=None, exclude_case=False):
        """
        Validate the object against a schema

        :param dict schema: a schema object to use to validate, it overrides the one
            that has been provided at object construction stage
        :param bool exclude_case: whether to exclude validated objects
            from the error. Useful when used with large configs
        """
        try:
            _validate(
                self.to_dict(expand=True), schema or getattr(self[IK], SCHEMA_KEY)
            )
        except ValidationError as e:
            _LOGGER.error(
                f"{self.__class__.__name__} object did not pass schema validation"
            )
            if getattr(self[IK], FILEPATH_KEY, None) is not None:
                # need to unlock locked files in case of validation error so that no
                # locks are left in place
                self.make_readonly()
            if not exclude_case:
                raise
            raise ValidationError(
                f"{self.__class__.__name__} object did not pass schema validation: "
                f"{e.message}"
            )
        _LOGGER.debug("Validated successfully")

    def write(self, filepath=None, schema=None, exclude_case=False):
        """
        Write the contents to a file.

        Make sure that the object has been created with write capabilities

        :param str filepath: a file path to write to
        :param dict schema: a schema object to use to validate, it overrides the one
            that has been provided at object construction stage
        :raise OSError: when the object has been created in a read only mode or other
            process has locked the file
        :raise TypeError: when the filepath cannot be determined. This takes place only
            if YacAttMap initialized with a Mapping as an input, not read from file.
        :raise OSError: when the write is called on an object with no write capabilities
            or when writing to a file that is locked by a different object
        :return str: the path to the created files
        """
        if getattr(self[IK], RO_KEY, False) and filepath is None:
            raise OSError(
                "You can't call write on an object that was created in read-only mode."
            )
        if schema is not None or getattr(self[IK], WRITE_VALIDATE_KEY):
            self.validate(schema=schema, exclude_case=exclude_case)
        filepath = _check_filepath(filepath or getattr(self[IK], FILEPATH_KEY, None))
        lock = make_lock_path(filepath)
        if filepath != getattr(self[IK], FILEPATH_KEY, None):
            if os.path.exists(filepath):
                if not os.path.exists(lock):
                    warnings.warn(
                        "Writing to a non-locked, existing file. Beware of collisions.",
                        UserWarning,
                    )
                else:
                    raise OSError(
                        f"The file '{filepath}' is locked by a different process"
                    )
            if getattr(self[IK], FILEPATH_KEY, None):
                self.make_readonly()
            setattr(self[IK], FILEPATH_KEY, filepath)
            create_lock(filepath, getattr(self[IK], WAIT_MAX_KEY, DEFAULT_WAIT_TIME))
        setattr(self[IK], RO_KEY, False)
        with open(filepath, "w") as f:
            f.write(self.to_yaml())
        abs_path = os.path.abspath(filepath)
        _LOGGER.debug(f"Wrote to a file: {abs_path}")
        return os.path.abspath(abs_path)

    @staticmethod
    def _remove_lock(filepath):
        """
        Remove lock

        :param str filepath: path to the file to remove the lock for. Not the
            path to the lock!
        :return bool: whether the lock was found and removed
        """
        lock = make_lock_path(_check_filepath(filepath))
        if os.path.exists(lock):
            os.remove(lock)
            return True
        return False

    def make_readonly(self):
        """
        Remove lock and make the object read only.

        :return bool: a logical indicating whether any locks were removed
        """
        if self._remove_lock(getattr(self[IK], FILEPATH_KEY, None)):
            setattr(self[IK], RO_KEY, True)
            _LOGGER.debug("Made object read-only")
            return True
        return False

    def make_writable(self, filepath=None):
        """
        Grant write capabilities to the object and re-read the file.

        Any changes made to the attributes are overwritten so that the object
        reflects the contents of the specified config file

        :param str filepath: path to the file that the contents will be written to
        :return YacAttMap: updated object
        """
        if not getattr(self[IK], RO_KEY, True):
            _LOGGER.info(
                "Object is already writable, path: {}".format(
                    getattr(self[IK], FILEPATH_KEY, None)
                )
            )
            return self
        ori_fp = getattr(self[IK], FILEPATH_KEY, None)
        if filepath and ori_fp != filepath:
            # file path has changed, unlock the previously used file if exists
            if ori_fp:
                self._remove_lock(ori_fp)
        filepath = _check_filepath(filepath or ori_fp)
        create_lock(filepath, getattr(self[IK], WAIT_MAX_KEY, DEFAULT_WAIT_TIME))
        try:
            self._reinit(filepath)
        except OSError:
            _LOGGER.debug("File '{}' not found".format(filepath))
            pass
        except Exception as e:
            self._reinit()
            _LOGGER.info(
                "File '{}' was not read, got an exception: {}".format(filepath, e)
            )
        setattr(self[IK], RO_KEY, False)
        setattr(self[IK], FILEPATH_KEY, filepath)
        _LOGGER.debug("Made object writable")
        return self

    @property
    def file_path(self):
        """
        Return the path to the config file or None if not set

        :return str | None: path to the file the object will would to
        """
        _warn_deprecated(obj=self)
        return getattr(self[IK], FILEPATH_KEY, None)

    @property
    def _file_path(self):
        """
        Return the path to the config file or None if not set

        :return str | None: path to the file the object will would to
        """
        _warn_deprecated(obj=self)
        return getattr(self[IK], FILEPATH_KEY, None)

    @property
    def writable(self):
        """
        Return writability flag or None if not set

        :return bool | None: whether the object is writable now
        """
        _warn_deprecated(obj=self)
        attr = getattr(self[IK], RO_KEY, None)
        return attr if attr is None else not attr


def _warn_deprecated(obj):
    fun_name = _getframe().f_back.f_code.co_name
    warnings.warn(
        f"The '{fun_name}' property is deprecated and will be removed in a future release."
        f' Use {obj.__class__.__name__}["{IK}"]["{fun_name}"] instead.',
        UserWarning,
        stacklevel=4,
    )


def _check_filepath(filepath):
    """
    Validate if the filepath attr/arg is a str

    :param str filepath: object to validate
    :return str: validated filepath
    :raise TypeError: if the filepath is not a string
    """
    # might be useful if we want to have multiple locked paths in the future
    # def _check_string(obj):
    #     """ check if object is a string or a list of strings """
    #     return bool(obj) and all(isinstance(elem, str) for elem in obj)
    if not isinstance(filepath, str):
        raise TypeError(
            f"No valid filepath provided. It has to be a str, got: {filepath.__class__.__name__}"
        )
    return filepath


def load_yaml(filepath: Union[str, Path]) -> dict:
    """
    Load a local or remote YAML file into a Python dict

    :param str filepath: path to the file to read
    :raises ConnectionError: if the remote YAML file reading fails
    :return dict: loaded yaml data
    """
    if is_url(filepath):
        _LOGGER.debug(f"Got URL: {filepath}")
        try:
            response = urlopen(filepath)
        except Exception as e:
            raise ConnectionError(
                f"Could not load remote file: {filepath}. "
                f"Original exception: {getattr(e, 'message', repr(e))}"
            )
        else:
            data = response.read().decode("utf-8")
            return yaml.safe_load(data)
    else:
        with open(os.path.abspath(filepath), "r") as f:
            data = yaml.safe_load(f)
        return data


def select_config(
    config_filepath: str = None,
    config_env_vars=None,
    default_config_filepath: str = None,
    check_exist: bool = True,
    on_missing=lambda fp: IOError(fp),
    strict_env: bool = False,
    config_name=None,
) -> str:
    """
    Selects the config file to load.

    This uses a priority ordering to first choose a config filepath if it's given,
    but if not, then look in a priority list of environment variables and choose
    the first available filepath to return.

    :param str | NoneType config_filepath: direct filepath specification
    :param Iterable[str] | NoneType config_env_vars: names of environment
        variables to try for config filepaths
    :param str default_config_filepath: default value if no other alternative
        resolution succeeds
    :param bool check_exist: whether to check for path existence as file
    :param function(str) -> object on_missing: what to do with a filepath if it
        doesn't exist
    :param bool strict_env: whether to raise an exception if no file path provided
        and environment variables do not point to any files
    raise: OSError: when strict environment variables validation is not passed
    """

    # First priority: given file
    if type(config_name) is str:
        config_name = f"{config_name} "
    else:
        config_name = ""

    if type(config_env_vars) is str:
        config_env_vars = [config_env_vars]

    if config_filepath:
        config_filepath = os.path.expandvars(config_filepath)
        if not check_exist or os.path.isfile(config_filepath):
            return os.path.abspath(config_filepath)
        _LOGGER.error(f"{config_name}config file path isn't a file: {config_filepath}")
        result = on_missing(config_filepath)
        if isinstance(result, Exception):
            raise result
        return os.path.abspath(result)

    _LOGGER.debug(f"No local {config_name}config file was provided.")
    selected_filepath = None

    # Second priority: environment variables (in order)
    if config_env_vars:
        _LOGGER.debug(
            f"Checking environment variables '{config_env_vars}' for {config_name}config"
        )

        for env_var in config_env_vars:
            result = os.environ.get(env_var, None)
            if result == None:
                _LOGGER.debug(f"Env var '{env_var}' not set.")
                continue
            elif result == "":
                _LOGGER.debug(f"Env var '{env_var}' exists, but value is empty.")
                continue
            elif not os.path.isfile(result):
                _LOGGER.debug(f"Env var '{env_var}' file not found: {result}")
                continue
            else:
                _LOGGER.debug(f"Found {config_name}config file in {env_var}: {result}")
                selected_filepath = result

    if selected_filepath is None:
        # Third priority: default filepath
        if default_config_filepath:
            _LOGGER.info(
                f"Using default {config_name}config. You may specify in env var: {str(config_env_vars)}"
            )
            return default_config_filepath
        else:
            if strict_env:
                raise OSError("Unable to select config file.")

            _LOGGER.info(f"Could not locate {config_name}config file.")
            return None
    return (
        os.path.abspath(selected_filepath) if selected_filepath else selected_filepath
    )


def deep_update(old, new):
    """
    Recursively update nested dict, modifying source
    """
    for key, value in new.items():
        if isinstance(value, Mapping) and value:
            old[key] = deep_update(old.get(key, {}), value)
        else:
            old[key] = new[key]
    return old
