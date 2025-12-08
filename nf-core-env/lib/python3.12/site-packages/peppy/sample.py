import glob
import os
from collections.abc import Mapping
from copy import copy as cp
from logging import getLogger
from string import Formatter
from typing import Optional, Union

import pandas as pd
import yaml
from pandas import Series, isnull

from .const import (
    CONFIG_FILE_KEY,
    DERIVED_SOURCES_KEY,
    PKG_NAME,
    PRJ_REF,
    SAMPLE_EDIT_FLAG_KEY,
    SAMPLE_NAME_ATTR,
    SAMPLE_SHEET_KEY,
)
from .exceptions import InvalidSampleTableFileException
from .utils import copy, grab_project_data
from .simple_attr_map import SimpleAttMap

_LOGGER = getLogger(PKG_NAME)


class SafeDict(dict):
    def __missing__(self, key):
        return "{" + key + "}"


@copy
class Sample(SimpleAttMap):
    """
    Class to model Samples based on a pandas Series.

    :param Mapping | pandas.core.series.Series series: Sample's data.
    """

    def __init__(self, series, prj=None):
        super(Sample, self).__init__()

        data = dict(series)
        _LOGGER.debug("Sample data: {data}")

        # Attach Project reference
        try:
            data_proj = data.pop(PRJ_REF)
        except (AttributeError, KeyError):
            data_proj = None

        self.update(**data)

        if data_proj and PRJ_REF not in self:
            self[PRJ_REF] = data_proj

        if PRJ_REF in self and prj:
            _LOGGER.warning(
                "Project data provided both in data and as separate"
                " constructor argument; using direct argument"
            )
        if prj:
            self[PRJ_REF] = prj
        if not self.get(PRJ_REF):
            # Force empty attmaps to null and ensure something's set.
            self[PRJ_REF] = None
            _LOGGER.debug("No project reference for sample")
        else:
            prefix = "Project reference on a sample must be an instance of dict"

            if not isinstance(self[PRJ_REF], Mapping):
                raise TypeError(f"{prefix}; got {type(self[PRJ_REF]).__name__}")
        self._derived_cols_done = []
        self._attributes = list(series.keys())

    def get_sheet_dict(self):
        """
        Create a K-V pairs for items originally passed in via the sample sheet.
        This is useful for summarizing; it provides a representation of the
        sample that excludes things like config files and derived entries.

        :return OrderedDict: mapping from name to value for data elements
            originally provided via the sample sheet (i.e., the a map-like
            representation of the instance, excluding derived items)
        """
        return dict([[k, self[k]] for k in self._attributes])

    def to_dict(self, add_prj_ref=False):
        """
        Serializes itself as dict object.

        :param bool add_prj_ref: whether the project reference bound do the
            Sample object should be included in the YAML representation
        :return dict: dict representation of this Sample
        """

        def _obj2dict(obj, name=None):
            """
            Build representation of object as a dict, recursively
            for all objects that might be attributes of self.

            :param object obj: what to serialize to write to YAML.
            :param str name: name of the object to represent.
            :param Iterable[str] to_skip: names of attributes to ignore.
            """
            if name:
                _LOGGER.log(5, "Converting to dict: {name}")
            if isinstance(obj, list):
                return [_obj2dict(i) for i in obj]
            elif isinstance(obj, Mapping):
                return {
                    k: _obj2dict(v, name=k)
                    for k, v in obj.items()
                    if not k.startswith("_")
                }
            if isinstance(obj, set):
                return [_obj2dict(i) for i in obj]
            elif isinstance(obj, Series):
                _LOGGER.warning("Serializing series as mapping, not array-like")
                return obj.to_dict()
            elif hasattr(obj, "dtype"):  # numpy data types
                # TODO: this fails with ValueError for multi-element array.
                return obj.item()
            elif isnull(obj):
                # Missing values as evaluated by pandas.isnull().
                # This gets correctly written into yaml.
                return None
            else:
                return obj

        serial = _obj2dict(self)
        if add_prj_ref:
            serial.update({"prj": grab_project_data(self[PRJ_REF])})
        return serial

    def to_yaml(
        self, path: Optional[str] = None, add_prj_ref=False
    ) -> Union[str, None]:
        """
        Serializes itself in YAML format. Writes to file if path is provided, else returns string representation.

        :param str path: A file path to write yaml to; provide this or
            the subs_folder_path, defaults to None
        :param bool add_prj_ref: whether the project reference bound do the
            Sample object should be included in the YAML representation
        :return str | None: returns string representation of sample yaml or None
        """
        serial = self.to_dict(add_prj_ref=add_prj_ref)
        if path:
            path = os.path.expandvars(path)
            if os.path.exists(os.path.dirname(path)):
                with open(path, "w") as outfile:
                    try:
                        yaml_data = yaml.safe_dump(serial, default_flow_style=False)
                    except yaml.representer.RepresenterError:
                        _LOGGER.error("Serialized sample data: {}".format(serial))
                        raise
                    outfile.write(yaml_data)
                    _LOGGER.debug("Sample data written to: {}".format(path))
            else:
                _LOGGER.warning(
                    "Could not write sample data to: {}. "
                    "Directory does not exist".format(path)
                )
                return
        else:
            yaml_data = yaml.safe_dump(serial, stream=None, default_flow_style=False)
            return yaml_data

    def derive_attribute(self, data_sources, attr_name):
        """
        Uses the template path provided in the project config section
        "data_sources" to piece together an actual path by substituting
        variables (encoded by "{variable}"") with sample attributes.

        :param Mapping data_sources: mapping from key name (as a value in
            a cell of a tabular data structure) to, e.g., filepath
        :param str attr_name: Name of sample attribute
            (equivalently, sample sheet column) specifying a derived column.
        :return str: regex expansion of data source specified in configuration,
            with variable substitutions made
        :raises ValueError: if argument to data_sources parameter is null/empty
        """

        def _format_regex(regex, items):
            """
            Format derived source with object attributes

            :param str regex: string to format,
                e.g. {identifier}{file_id}_data.txt
            :param Iterable[Iterable[Iterable | str]] items: items to format
                the string with
            :raise InvalidSampleTableFileException: if after merging
                subannotations the lengths of multi-value attrs are not even
            :return Iterable | str: formatted regex string(s)
            """
            keys = [i[1] for i in Formatter().parse(regex) if i[1] is not None]
            if not keys:
                return [regex]
            if "$" in regex:
                _LOGGER.warning(
                    "Not all environment variables were populated "
                    "in derived attribute source: {}".format(regex)
                )
            attr_lens = [
                len(v) for k, v in items.items() if (isinstance(v, list) and k in keys)
            ]
            if not bool(attr_lens):
                return [_safe_format(regex, items)]
            if len(set(attr_lens)) != 1:
                msg = (
                    "All attributes to format the {} ({}) have to be the "
                    "same length, got: {}. Correct your {}".format(
                        DERIVED_SOURCES_KEY, regex, attr_lens, SAMPLE_SHEET_KEY
                    )
                )
                raise InvalidSampleTableFileException(msg)
            vals = []
            for i in range(0, attr_lens[0]):
                items_cpy = cp(items)
                for k in keys:
                    if isinstance(items_cpy[k], list):
                        items_cpy[k] = items_cpy[k][i]
                vals.append(_safe_format(regex, items_cpy))
            return vals

        def _safe_format(s, values):
            """
            Safely format string.

            If the values are missing the key is wrapped in curly braces.
            This is intended to preserve the environment variables specified
            using curly braces notation, for example: "${ENVVAR}/{sample_attr}"
            would result in "${ENVVAR}/populated" rather than a KeyError.

            :param str s: string with curly braces placeholders to populate
            :param Mapping values: key-value pairs to pupulate string with
            :return str: populated string
            """
            return Formatter().vformat(s, (), SafeDict(values))

        def _glob_regex(patterns):
            """
            Perform unix style pathname pattern expansion for multiple patterns

            :param Iterable[str] patterns: patterns to expand
            :return str | Iterable[str]: expanded patterns
            """
            outputs = []
            for p in patterns:
                if "*" in p or "[" in p:
                    _LOGGER.debug("Pre-glob: {}".format(p))
                    val_globbed = sorted(glob.glob(p))
                    if not val_globbed:
                        _LOGGER.debug("No files match the glob: '{}'".format(p))
                    else:
                        p = val_globbed
                        _LOGGER.debug("Post-glob: {}".format(p))

                outputs.extend(p if isinstance(p, list) else [p])
            return outputs if len(outputs) > 1 else outputs[0]

        if not data_sources:
            return None
        sn = self[SAMPLE_NAME_ATTR] if SAMPLE_NAME_ATTR in self else "this sample"
        try:
            source_key = self[attr_name]
        except AttributeError:
            reason = (
                "'{attr}': to locate sample's derived attribute source, "
                "provide the name of a key from '{sources}' or ensure "
                "sample has attribute '{attr}'".format(
                    attr=attr_name, sources=DERIVED_SOURCES_KEY
                )
            )
            raise AttributeError(reason)

        try:
            regex = data_sources[source_key]
            _LOGGER.debug("Data sources: {data_sources}")
        except KeyError:
            _LOGGER.debug(
                f"{sn}: config lacks entry for {DERIVED_SOURCES_KEY} key: "
                f"'{source_key}' in column '{attr_name}'; known: {data_sources.keys()}"
            )
            return ""
        deriv_exc_base = (
            f"In sample '{sn}' cannot correctly parse derived "
            f"attribute source: {regex}."
        )
        try:
            expanded_regex = os.path.expandvars(regex)
            vals = _format_regex(expanded_regex, dict(self.items()))
            _LOGGER.debug("Formatted regex: {}".format(vals))
        except KeyError as ke:
            _LOGGER.warning(f"{deriv_exc_base} Can't access {str(ke)} attribute")
        except Exception as e:
            _LOGGER.warning(f"{deriv_exc_base} Caught exception: {str(e)}")
        else:
            return _glob_regex(vals)
        return None

    @property
    def project(self):
        """
        Get the project mapping

        :return peppy.Project: project object the sample was created from
        """
        return self[PRJ_REF]

    @property
    def sample_name(self):
        """
        Get the sample's name

        :return str: current sample name derived from project's st_index
        """

        return self[self[PRJ_REF].st_index]

    # The __reduce__ function provides an interface for
    # correct object serialization with the pickle module.
    def __reduce__(self):
        return (
            self.__class__,
            (pd.Series(self.to_dict()),),
            # (None, {}),
            # iter([]),
            # iter({PRJ_REF: self[PRJ_REF]}.items()),
        )

    def __len__(self):
        return len(self.to_dict())

    def __str__(self, max_attr=10):
        """Representation in interpreter."""
        if len(self) == 0:
            return ""
        head = "Sample"
        try:
            head += f" '{self[SAMPLE_NAME_ATTR]}'"
        except KeyError:
            pass
        try:
            prj_cfg = self[PRJ_REF][CONFIG_FILE_KEY]
        except (KeyError, TypeError):
            pass
        else:
            head += f" in Project ({prj_cfg})"
        pub_attrs = {k: v for k, v in self.items() if not k.startswith("_")}
        maxlen = max(map(len, pub_attrs.keys())) + 2
        attrs = ""
        counter = 0
        for k, v in pub_attrs.items():
            key_to_show = (k + ":").ljust(maxlen)
            if not isinstance(v, list):
                val_to_show = v
            else:
                try:
                    val_to_show = ", ".join([i for i in v if v is not None])
                except TypeError:
                    val_to_show = "None"
            attrs += f"\n{key_to_show}{val_to_show}"
            if counter == max_attr:
                attrs += "\n\n...".ljust(maxlen) + f"(showing first {max_attr})"
                break
            counter += 1
        return head + "\n" + attrs

    def __repr__(self):
        return str(self)

    def _excl_from_eq(self, k):
        """Exclude the Project reference from object comparison."""
        return k == PRJ_REF or super(Sample, self)._excl_from_eq(k)

    def _excl_from_repr(self, k, cls):
        """Exclude the Project reference from representation."""
        return k.startswith("_") or super(Sample, self)._excl_from_repr(k, cls)

    def _excl_classes_from_todict(self):
        """Exclude pandas.DataFrame from dict representation"""
        return (pd.DataFrame,)

    def _try_touch_samples(self):
        """
        Safely sets sample edited flag to true
        """
        try:
            self[PRJ_REF][SAMPLE_EDIT_FLAG_KEY] = True
        except (KeyError, AttributeError, TypeError):
            pass
