import logging
import os
import yaml

from collections.abc import Iterable, Mapping
from jsonschema import validate as _validate
from jsonschema.exceptions import ValidationError
from sys import _getframe
from ubiquerg import (
    expandpath,
    is_url,
    ThreeLocker,
    ensure_locked,
    locked_read_file,
    READ,
    WRITE,
)

from ._version import __version__

_LOGGER = logging.getLogger(__name__)
_LOGGER.debug(f"Using yacman version {__version__}")

# Hack for yaml string indexes
# Credit: Anthon
# https://stackoverflow.com/questions/50045617
# https://stackoverflow.com/questions/5121931
# The idea is: if you have yaml keys that can be interpreted as an int or a float,
# then the yaml loader will convert them into an int or a float, and you would
# need to access them with dict[2] instead of dict['2']. But since we always
# expect the keys to be strings, this doesn't work. So, here we are adjusting
# the loader to keep everything as a string.

# Only do once.
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

# Constants: to do, remove these

DEFAULT_WAIT_TIME = 60
# LOCK_PREFIX = "lock."
SCHEMA_KEY = "schema"
FILEPATH_KEY = "file_path"

from collections.abc import MutableMapping

# Since read and write are now different context managers, we have to
# separate them like this, instead of using __enter__ and __exit__ on the class
# itself, which only allows one type of context manager.


class FutureYAMLConfigManager(MutableMapping):
    """
    A YAML configuration manager, providing file locking, loading,
    writing, etc.  for YAML configuration files.
    """

    def __init__(
        self,
        entries=None,
        wait_max=DEFAULT_WAIT_TIME,
        strict_ro_locks=False,
        schema_source=None,
        validate_on_write=False,
    ):
        """
        Object constructor

        :param Iterable[(str, object)] | Mapping[str, object] entries: YAML collection
            of key-value pairs.
        :param str yamldata: YAML-formatted string
        :param int wait_max: how long to wait for creating an object when the file
            that data will be read from is locked
        :param bool strict_ro_locks: By default, we allow RO filesystems that can't be locked.
            Turn on strict_ro_locks to error if locks cannot be enforced on readonly filesystems.
        :param bool skip_read_lock: whether the file should not be locked for reading
            when object is created in read only mode
        :param str schema_source: path or a URL to a jsonschema in YAML format to use
            for optional config validation. If this argument is provided the object
            is always validated at least once, at the object creation stage.
        :param bool validate_on_write: a boolean indicating whether the object should be
            validated every time the `write` method is executed, which is
            a way of preventing invalid config writing

        """

        # Settings for this config object
        self.filepath = None
        self.wait_max = wait_max
        self.schema_source = schema_source
        self.validate_on_write = validate_on_write
        self.strict_ro_locks = strict_ro_locks
        self.locker = None

        # We store the values in a dict under .data
        if isinstance(entries, list):
            self.data = entries
        else:
            self.data = dict(entries or {})
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
            setattr(self, SCHEMA_KEY, load_yaml(sp))
            self.validate()

    @classmethod
    def from_obj(cls, entries: object, **kwargs):
        """
        Initialize from a Python object (dict, list, or primitive).

        :param obj entries: object to initialize from.
        :param kwargs: Keyword arguments to pass to the constructor.
        """
        return cls(entries, **kwargs)

    @classmethod
    def from_yaml_data(cls, yamldata, **kwargs):
        """
        Initialize from a YAML string.

        :param str yamldata: YAML-formatted string.
        :param kwargs: Keyword arguments to pass to the constructor.
        """
        entries = yaml.load(yamldata, yaml.SafeLoader)
        return cls(entries, **kwargs)

    @classmethod
    def from_yaml_file(cls, filepath: str, create_file: bool = False, **kwargs):
        """
        Initialize from a YAML file.

        :param str filepath: Path to the YAML config file.
        :param str create_file: Create a file at filepath if it doesn't exist.
        :param kwargs: Keyword arguments to pass to the constructor.
        """

        file_contents = locked_read_file(filepath, create_file=create_file)
        entries = yaml.load(file_contents, yaml.SafeLoader)
        ref = cls(entries, **kwargs)
        ref.locker = ThreeLocker(filepath)
        ref.filepath = filepath
        return ref

    def update_from_yaml_file(self, filepath=None):
        if filepath is not None:  # set filepath to update filepath if uninitialized
            if self.filepath is not None:
                self.filepath = filepath
        self.data.update(load_yaml(filepath))
        return

    def update_from_yaml_data(self, yamldata=None):
        self.data.update(yaml.load(yamldata, yaml.SafeLoader))
        return

    def update_from_obj(self, entries=None):
        self.data.update(entries)
        return

    @property
    def settings(self):
        return {
            "wait_max": self.wait_max,
            "schema_source": self.schema_source,
            "validate_on_write": self.validate_on_write,
            "locked": self.locked,
            "strict_ro_locks": self.strict_ro_locks,
        }

    def __del__(self):
        if hasattr(self, "locker"):
            del self.locker

    def __repr__(self):
        # Render the data in a nice way
        return self.to_yaml(self.data)

    def __enter__(self):
        raise NotImplementedError(
            "Use the 'read_lock' and 'write_lock' context managers."
        )

    def __exit__(self):
        raise NotImplementedError(
            "Use the 'read_lock' and 'write_lock' context managers."
        )

    @ensure_locked(READ)
    def rebase(self, filepath=None):
        """
        Reload the object from file, then update with current information

        :param str filepath: path to the file that should be read
        """
        fp = filepath or self.locker.filepath
        if fp is not None:
            local_data = self.data
            self.data = load_yaml(fp)
            _LOGGER.debug(f"Rebased {local_data} with {self.data} from {fp}")
            if self.data is None:
                self.data = local_data
            else:
                deep_update(self.data, local_data)
        else:
            _LOGGER.warning("Rebase has no effect if no filepath given")

        return self

    @ensure_locked(READ)
    def reset(self, filepath=None):
        """
        Reset dict contents to file contents, or to empty dict if no filepath found.
        """
        fp = filepath or self.locker.filepath
        if fp is not None:
            self.data = load_yaml(fp)
        else:
            self.data = {}
        return self

    def validate(self, schema=None, exclude_case=False):
        """
        Validate the object against a schema

        :param dict schema: a schema object to use to validate, it overrides the one
            that has been provided at object construction stage
        :param bool exclude_case: whether to exclude validated objects
            from the error. Useful when used with large configs
        """
        try:
            _validate(self.to_dict(expand=True), schema or getattr(self, SCHEMA_KEY))
        except ValidationError as e:
            _LOGGER.error(
                f"{self.__class__.__name__} object did not pass schema validation"
            )
            # if getattr(self, FILEPATH_KEY, None) is not None:
            # need to unlock locked files in case of validation error so that no
            # locks are left in place
            # self.make_readonly()
            # commented out because I think this is taken care of my context managers now
            if not exclude_case:
                raise
            raise ValidationError(
                f"{self.__class__.__name__} object did not pass schema validation: "
                f"{e.message}"
            )
        _LOGGER.debug("Validated successfully")

    @ensure_locked(WRITE)
    def write(self, schema=None, exclude_case=False):
        """
        Write the contents to the file backing this object.

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
        if not self.locker.filepath:
            raise OSError("Must provide a filepath to write.")

        _check_filepath(self.locker.filepath)
        _LOGGER.debug(f"Writing to file '{self.locker.filepath}'")
        with open(self.locker.filepath, "w") as f:
            f.write(self.to_yaml())

        if schema is not None or self.validate_on_write:
            self.validate(schema=schema, exclude_case=exclude_case)

        abs_path = os.path.abspath(self.locker.filepath)
        _LOGGER.debug(f"Wrote to a file: {abs_path}")
        return os.path.abspath(abs_path)

    def write_copy(self, filepath=None):
        """
        Write the contents to an external file.

        :param str filepath: a file path to write to
        """

        _LOGGER.debug(f"Writing to file '{filepath}'")
        with open(filepath, "w") as f:
            f.write(self.to_yaml())
        return filepath

    def to_yaml(self, trailing_newline=False, expand=False):
        """
        Get text for YAML representation.

        :param bool trailing_newline: whether to add trailing newline
        :param bool expand: whether to expand paths in values
        :return str: YAML text representation of this instance.
        """

        if expand:
            return yaml.dump(self.exp, default_flow_style=False)
        return yaml.dump(self.data, default_flow_style=False) + (
            "\n" if trailing_newline else ""
        )

    def to_dict(self, expand=True):
        # Seems like it's probably not necessary; can just use the object now.
        # but for backwards compatibility.
        return self.data

    def __setitem__(self, item, value):
        self.data[item] = value

    def __getitem__(self, item):
        """
        Fetch the value of given key.

        :param hashable item: key for which to fetch value
        :return object: value mapped to given key, if available
        :raise KeyError: if the requested key is unmapped.
        """
        return self.data[item]

    @property
    def exp(self) -> dict:
        """
        Returns a copy of the object's data elements with env vars and user vars
        expanded. Use it like: object.exp["item"]
        """
        return _safely_expand_path(self.data)

    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)

    def __delitem__(self, key):
        value = self[key]
        del self.data[key]
        self.pop(value, None)

    def priority_get(
        self,
        arg_name: str,
        env_var: str = None,
        default: str = None,
        override: str = None,
        strict: bool = False,
    ):
        """
        Helper function to select a value from a config, or, if missing, then go to an env var.

        :param str arg_name: Argument to retrieve
        :param bool strict: Should missing args raise an error? False=warning
        :param env_var: Env var to retrieve from should it be missing from the cfg
        """
        if override:
            return override
        if self.data.get(arg_name) is not None:
            return self.data[arg_name]
        if env_var is not None:
            arg = os.getenv(env_var, None)
            if arg is not None:
                _LOGGER.debug(f"Value '{arg}' sourced from '{env_var}' env var")
                return expandpath(arg)
        if default is not None:
            return default
        if strict:
            message = (
                "Value for required argument '{arg_name}' could not be determined."
            )
            _LOGGER.warning(message)
            raise Exception(message)


# A big issue here is: if you route the __getitem__ through this,
# then it returns a copy of the data, rather than the data itself.
# That's the point, so we don't adjust it. But then you can't use multi-level
# item setting, like ycm["x"]["y"] = value, because ycm['x'] returns a different
# dict, and so you're updating that copy of it.
# The solution is that we have to route expansion through a separate property,
# so the setitem syntax can remain intact while preserving original values.
def _safely_expand_path(x):
    if isinstance(x, str):
        return expandpath(x)
    elif isinstance(x, Mapping):
        return {k: _safely_expand_path(v) for k, v in x.items()}
    return x


def _unsafely_expand_path(x):
    if isinstance(x, str):
        return expandpath(x)
    elif isinstance(x, Mapping):
        for k in x.keys():
            x[k] = _safely_expand_path(x[k])
        return x
        # return {k: _safely_expand_path(v) for k, v in x.items()}
    return x


def _check_filepath(filepath):
    """
    Validate if the filepath is a str

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
            f"No valid filepath provided. It must be a str, got: {filepath.__class__.__name__}"
        )
    return filepath


def _warn_deprecated(obj):
    fun_name = _getframe().f_back.f_code.co_name
    warnings.warn(
        f"The '{fun_name}' property is deprecated and will be removed in a future release."
        f' Use {obj.__class__.__name__}["{IK}"]["{fun_name}"] instead.',
        UserWarning,
        stacklevel=4,
    )


def load_yaml(filepath):
    """Load a yaml file into a python dict"""

    def read_yaml_file(filepath):
        """
        Read a YAML file

        :param str filepath: path to the file to read
        :return dict: read data
        """
        with open(filepath, "r") as f:
            data = yaml.safe_load(f)
        return data

    if is_url(filepath):
        _LOGGER.debug(f"Got URL: {filepath}")
        try:  # python3
            from urllib.error import HTTPError
            from urllib.request import urlopen
        except:  # python2
            from urllib2 import URLError as HTTPError
            from urllib2 import urlopen
        try:
            response = urlopen(filepath)
        except HTTPError as e:
            raise e
        data = response.read()  # a `bytes` object
        text = data.decode("utf-8")
        return yaml.safe_load(text)
    else:
        return read_yaml_file(filepath)


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
