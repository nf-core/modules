"""Abstraction of a parse of a schema definition"""

import copy
import logging
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Union

import yacman

from .const import (
    CANONICAL_TYPES,
    CLASSES_BY_TYPE,
    RECORD_IDENTIFIER,
    SCHEMA_DESC_KEY,
    SCHEMA_ITEMS_KEY,
    SCHEMA_PROP_KEY,
    SCHEMA_TYPE_KEY,
)
from .exceptions import SchemaError

_LOGGER = logging.getLogger(__name__)


NULL_MAPPING_VALUE = {}
SCHEMA_PIPELINE_NAME_KEY = "pipeline_name"


# The columns associated with the file and image types
PATH_COL_SPEC = (Path, ...)


def _safe_pop_one_mapping(
    mappingkey: str,
    data: Dict[str, Any],
    info_name: str,
    subkeys: Optional[List[str]] = None,
) -> Any:
    """
    mapping key: the dict key where the sample, project or status values are stored, e.g. data["mappingkey"]
    subkeys: if using JSON schema, the dict is nested further, e.g. data["properties"]["samples"]["mappingkey"]
    """
    if subkeys:
        try:
            value = data[subkeys[0]].pop(mappingkey, NULL_MAPPING_VALUE)
        except KeyError:
            value = {}
        if isinstance(value, Mapping):
            return value
    else:
        value = data.pop(mappingkey, NULL_MAPPING_VALUE)
        if isinstance(value, Mapping):
            return value
        raise SchemaError(
            f"{info_name} info in schema definition has invalid type: {type(value).__name__}"
        )


class ParsedSchema(object):
    """
    Store the results of parsing a pipestat schema configuration file.

    In particular, there are different 'levels' (concepts, really) at which schema
    elements may be defined; namely, there may be project-, sample-, or status-related
    schema information in a configuration file.

    This class tames this complexity relative to interacting directly with a raw
    Mapping-like object that would result from a parse, providing accessors for each
    of the key groupings of schema information, as well as the name of the pipeline
    for which the schema is written.
    """

    _PROJECT_KEY = "project"
    _SAMPLES_KEY = "samples"
    _STATUS_KEY = "status"

    def __init__(self, data: Union[Dict[str, Any], Path, str]) -> None:
        # initial validation and parse
        if not isinstance(data, dict):
            data = yacman.load_yaml(data)

        # Keep a copy of the original schema
        self.original_schema = copy.deepcopy(data)
        # Create a copy of a resolved schema
        if "$defs" in data.keys():
            self.resolved_schema = replace_JSON_refs(copy.deepcopy(data), data)
            self.resolved_schema.pop("$defs")
        else:
            self.resolved_schema = copy.deepcopy(data)

        data = copy.deepcopy(data)

        # Currently supporting backwards compatibility with old output schema while now also supporting a JSON schema:
        if "properties" in list(data.keys()):
            # Assume top-level properties key implies proper JSON schema.

            self._pipeline_name = data["properties"].pop(SCHEMA_PIPELINE_NAME_KEY, None)

            # Two passes for sample-data as it is now nested under items per #204
            sample_data = _safe_pop_one_mapping(
                subkeys=["samples"],
                data=data["properties"],
                info_name="sample-level",
                mappingkey="items",
            )
            sample_data = _safe_pop_one_mapping(
                data=sample_data,
                info_name="sample-level",
                mappingkey="properties",
            )

            prj_data = _safe_pop_one_mapping(
                subkeys=["project"],
                data=data["properties"],
                info_name="project-level",
                mappingkey="properties",
            )

            self._status_data = _safe_pop_one_mapping(
                subkeys=["status"],
                data=data["properties"],
                info_name="status",
                mappingkey="properties",
            )
            # TODO We should add the ability to look at an external source beyond the source schema (data)
            self._sample_level_data = replace_JSON_refs(sample_data, data)

            self._project_level_data = replace_JSON_refs(prj_data, data)

        else:
            self._pipeline_name = data.pop(SCHEMA_PIPELINE_NAME_KEY, None)
            sample_data = _safe_pop_one_mapping(
                mappingkey=self._SAMPLES_KEY, data=data, info_name="sample-level"
            )
            prj_data = _safe_pop_one_mapping(
                mappingkey=self._PROJECT_KEY, data=data, info_name="project-level"
            )
            # Parse custom status declaration if present.
            self._status_data = _safe_pop_one_mapping(
                mappingkey=self._STATUS_KEY, data=data, info_name="status"
            )
            self._sample_level_data = _recursively_replace_custom_types(sample_data)

            self._project_level_data = _recursively_replace_custom_types(prj_data)

        if not isinstance(self._pipeline_name, str):
            raise SchemaError(
                f"Could not find valid pipeline identifier (key '{SCHEMA_PIPELINE_NAME_KEY}') in given schema data"
            )

        # Sample- and/or project-level data must be declared.
        if not self._sample_level_data and not self._project_level_data:
            raise SchemaError("Neither sample-level nor project-level data items are declared.")

        # Check that no reserved keywords were used as data items.
        resv_kwds = {"id", RECORD_IDENTIFIER}
        reserved_keywords_used = set()
        for data in [self.project_level_data, self.sample_level_data, self.status_data]:
            reserved_keywords_used |= set(data.keys()) & resv_kwds
        if reserved_keywords_used:
            raise SchemaError(
                f"{len(reserved_keywords_used)} reserved keyword(s) used: {', '.join(reserved_keywords_used)}"
            )

        # Check that no data item name overlap exists between project- and sample-level data.
        project_sample_overlap = set(self.project_level_data) & set(self.sample_level_data)
        if project_sample_overlap:
            raise SchemaError(
                f"Overlap between project- and sample-level keys: {', '.join(project_sample_overlap)}"
            )

    def __str__(self):
        """
        Generate string representation of the object.

        :return str: string representation of the object
        """
        res = f"{self.__class__.__name__} ({self._pipeline_name})"

        def add_props(props):
            res = ""
            if len(props) == 0:
                res += "\n - None"
            else:
                for k, v in props:
                    res += f"\n - {k} : {v}"
            return res

        if self._project_level_data is not None:
            res += "\n Project-level properties:"
            res += add_props(self._project_level_data.items())
        if self._sample_level_data is not None:
            res += "\n Sample-level properties:"
            res += add_props(self._sample_level_data.items())
        if self._status_data is not None:
            res += "\n Status properties:"
            res += add_props(self._status_data.items())
        return res

    def __repr__(self):
        """
        Generate string representation of the object.

        :return str: string representation of the object
        """
        return self.__str__()

    @property
    def pipeline_name(self):
        """Return the declared name for the pipeline for which this schema's written."""
        return self._pipeline_name

    @property
    def project_level_data(self):
        """Return information relevant for a project-level pipeline."""
        return copy.deepcopy(self._project_level_data)

    @property
    def results_data(self):
        """Return union of sample- and project-level information."""
        return {**self.project_level_data, **self.sample_level_data}

    @property
    def sample_level_data(self):
        """Return information relevant for a sample-level pipeline."""
        return copy.deepcopy(self._sample_level_data)

    @property
    def status_data(self):
        """Return information relevant to pipeline status."""
        return copy.deepcopy(self._status_data)

    @property
    def project_table_name(self):
        """Return the name of the database table for project-level information."""
        return self._table_name("project")

    @property
    def sample_table_name(self):
        """Return the name of the database table for sample-level information."""
        return self._table_name("sample")

    @staticmethod
    def _get_data_type(type_name):
        t = CLASSES_BY_TYPE[type_name]
        # return ARRAY if t == list else t
        return t

    @property
    def file_like_table_name(self):
        return self._table_name("files")

    def to_dict(self) -> Dict[str, Any]:
        """Create simple dictionary representation of this instance."""
        data = {SCHEMA_PIPELINE_NAME_KEY: self.pipeline_name}
        for key, values in [
            (self._PROJECT_KEY, self.project_level_data),
            (self._SAMPLES_KEY, self.sample_level_data),
            (self._STATUS_KEY, self.status_data),
        ]:
            if values:
                data[key] = values
        return data

    def _table_name(self, suffix: str) -> str:
        return f"{self.pipeline_name}__{suffix}"


def _recursively_replace_custom_types(s: Dict[str, Any]) -> Dict[str, Any]:
    """
    Replace the custom types in pipestat schema with canonical types

    :param dict s: schema to replace types in
    :return dict: schema with types replaced
    """
    for k, v in s.items():
        missing_req_keys = [req for req in [SCHEMA_TYPE_KEY, SCHEMA_DESC_KEY] if req not in v]
        if missing_req_keys:
            raise SchemaError(
                f"Result '{k}' is missing required key(s): {', '.join(missing_req_keys)}"
            )
        curr_type_name = v[SCHEMA_TYPE_KEY]
        if curr_type_name == "object" and SCHEMA_PROP_KEY in s[k]:
            _recursively_replace_custom_types(s[k][SCHEMA_PROP_KEY])
        if curr_type_name == "array" and SCHEMA_ITEMS_KEY in s[k]:
            _recursively_replace_custom_types(s[k][SCHEMA_ITEMS_KEY][SCHEMA_PROP_KEY])
        try:
            curr_type_spec = CANONICAL_TYPES[curr_type_name]
        except KeyError:
            continue
        spec = s.setdefault(k, {})
        spec.setdefault(SCHEMA_PROP_KEY, {}).update(curr_type_spec[SCHEMA_PROP_KEY])
        spec.setdefault("required", []).extend(curr_type_spec["required"])
        spec[SCHEMA_TYPE_KEY] = curr_type_spec[SCHEMA_TYPE_KEY]
    return s


def replace_JSON_refs(
    target_schema: Dict[str, Any], source_schema: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Recursively search and replace the $refs if they exist in schema, target_schema, and if their corresponding $defs
    exist in source schema, source_schema. If $defs  exist in the target schema and target_schema is the same as
    source_schema then deepcopy should be used such that target_schema = copy.deepcopy(source_schema)

    :param dict target_schema: schema to replace types in
    :param dict source_schema: source schema
    :return dict target_schema: schema with types replaced
    """

    for k, v in list(target_schema.items()):
        if isinstance(v, dict):
            replace_JSON_refs(target_schema[k], source_schema)
        if "$ref" == k:
            split_value = v.split("/")
            if len(split_value) != 3:
                raise SchemaError(
                    msg=f"$ref exists in source schema but path,{v} ,not valid, e.g. '#/$defs/file' "
                )
            if split_value[1] in source_schema and split_value[2] in source_schema[split_value[1]]:
                result = source_schema[split_value[1]][split_value[2]]
            else:
                result = None
            if result is not None:
                for key, value in result.items():
                    target_schema.update({key: value})
                del target_schema["$ref"]
            else:
                raise SchemaError(
                    msg=f"Could not find {split_value[1]} and {split_value[2]} in $def"
                )
    return target_schema
