"""Helpers without an obvious logical home."""

import logging
import os
from typing import Dict, Mapping, Type, Union
from urllib.request import urlopen

import yaml
from ubiquerg import expandpath, is_url

from .const import CONFIG_KEY, SAMPLE_TABLE_INDEX_KEY, SUBSAMPLE_TABLE_INDEX_KEY
from .exceptions import RemoteYAMLError

_LOGGER = logging.getLogger(__name__)


def copy(obj):
    def copy(self):
        """
        Copy self to a new object.
        """
        from copy import deepcopy

        return deepcopy(self)

    obj.copy = copy
    return obj


def make_abs_via_cfg(maybe_relpath, cfg_path, check_exists=False):
    """Ensure that a possibly relative path is absolute."""
    if not isinstance(maybe_relpath, str):
        raise TypeError(
            "Attempting to ensure non-text value is absolute path: {} ({})".format(
                maybe_relpath, type(maybe_relpath)
            )
        )
    if os.path.isabs(maybe_relpath) or is_url(maybe_relpath):
        _LOGGER.debug("Already absolute")
        return maybe_relpath
    # Maybe we have env vars that make the path absolute?
    expanded = expandpath(maybe_relpath)
    if os.path.isabs(expanded):
        _LOGGER.debug("Expanded: {}".format(expanded))
        return expanded
    # Set path to an absolute path, relative to project config.
    config_dirpath = os.path.dirname(cfg_path)
    _LOGGER.debug("config_dirpath: {}".format(config_dirpath))
    abs_path = os.path.join(config_dirpath, maybe_relpath)
    _LOGGER.debug("Expanded and/or made absolute: {}".format(abs_path))
    if check_exists and not os.path.exists(abs_path):
        raise OSError(f"Path made absolute does not exist: {abs_path}")
    return abs_path


def grab_project_data(prj):
    """
    From the given Project, grab Sample-independent data.

    There are some aspects of a Project of which it's beneficial for a Sample
    to be aware, particularly for post-hoc analysis. Since Sample objects
    within a Project are mutually independent, though, each doesn't need to
    know about any of the others. A Project manages its, Sample instances,
    so for each Sample knowledge of Project data is limited. This method
    facilitates adoption of that conceptual model.

    :param Project prj: Project from which to grab data
    :return Mapping: Sample-independent data sections from given Project
    """
    if not prj:
        return {}

    try:
        return dict(prj[CONFIG_KEY])
    except KeyError:
        raise KeyError("Project lacks section '{}'".format(CONFIG_KEY))


def make_list(arg: Union[list, str], obj_class: Type) -> list:
    """
    Convert an object of predefined class to a list of objects of that class or
    ensure a list is a list of objects of that class

    :param list[obj] | obj arg: string or a list of strings to listify
    :param str obj_class: name of the class of intrest

    :return list: list of objects of the predefined class

    :raise TypeError: if a faulty argument was provided
    """

    def _raise_faulty_arg():
        raise TypeError(
            "Provided argument has to be a List[{o}] or a {o}, "
            "got '{a}'".format(o=obj_class.__name__, a=arg.__class__.__name__)
        )

    if isinstance(arg, obj_class):
        return [arg]
    elif isinstance(arg, list):
        if not all(isinstance(i, obj_class) for i in arg):
            _raise_faulty_arg()
        else:
            return arg
    else:
        _raise_faulty_arg()


def _expandpath(path: str):
    """
    Expand a filesystem path that may or may not contain user/env vars.

    :param str path: path to expand
    :return str: expanded version of input path
    """
    return os.path.expandvars(os.path.expanduser(path))


def expand_paths(x: dict) -> dict:
    """
    Recursively expand paths in a dict.

    :param dict x: dict to expand
    :return dict: dict with expanded paths
    """
    if isinstance(x, str):
        return expandpath(x)
    elif isinstance(x, Mapping):
        return {k: expand_paths(v) for k, v in x.items()}
    return x


def load_yaml(filepath):
    """
    Load a local or remote YAML file into a Python dict

    :param str filepath: path to the file to read
    :raises RemoteYAMLError: if the remote YAML file reading fails
    :return dict: read data
    """
    if is_url(filepath):
        _LOGGER.debug(f"Got URL: {filepath}")
        try:
            response = urlopen(filepath)
        except Exception as e:
            raise RemoteYAMLError(
                f"Could not load remote file: {filepath}. "
                f"Original exception: {getattr(e, 'message', repr(e))}"
            )
        else:
            data = response.read().decode("utf-8")
            return expand_paths(yaml.safe_load(data))
    else:
        with open(os.path.abspath(filepath), "r") as f:
            data = yaml.safe_load(f)
        return expand_paths(data)


def is_cfg_or_anno(file_path, formats=None):
    """
    Determine if the input file seems to be a project config file (based on the file extension).
    :param str file_path: file path to examine
    :param dict formats: formats dict to use. Must include 'config' and 'annotation' keys.
    :raise ValueError: if the file seems to be neither a config nor an annotation
    :return bool: True if the file is a config, False if the file is an annotation
    """
    formats_dict = formats or {
        "config": (".yaml", ".yml"),
        "annotation": (".csv", ".tsv"),
    }
    if file_path is None:
        return None
    if file_path.lower().endswith(formats_dict["config"]):
        _LOGGER.debug(f"Creating a Project from a YAML file: {file_path}")
        return True
    elif file_path.lower().endswith(formats_dict["annotation"]):
        _LOGGER.debug(f"Creating a Project from a CSV file: {file_path}")
        return False
    raise ValueError(
        f"File path '{file_path}' does not point to an annotation or config. "
        f"Accepted extensions: {formats_dict}"
    )


def extract_custom_index_for_sample_table(pep_dictionary: Dict):
    """Extracts a custom index for the sample table if it exists"""
    return (
        pep_dictionary[SAMPLE_TABLE_INDEX_KEY]
        if SAMPLE_TABLE_INDEX_KEY in pep_dictionary
        else None
    )


def extract_custom_index_for_subsample_table(pep_dictionary: Dict):
    """Extracts a custom index for the subsample table if it exists"""
    return (
        pep_dictionary[SUBSAMPLE_TABLE_INDEX_KEY]
        if SUBSAMPLE_TABLE_INDEX_KEY in pep_dictionary
        else None
    )
