"""Abstraction of a parse of a schema definition"""

import copy
import datetime
import logging
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

from pydantic import ConfigDict, create_model
from sqlalchemy import Column, null
from sqlalchemy.dialects.postgresql import JSONB
from sqlmodel import Field, SQLModel

from pipestat.const import (
    CANONICAL_TYPES,
    CLASSES_BY_TYPE,
    CREATED_TIME,
    ID_KEY,
    MODIFIED_TIME,
    PIPELINE_NAME,
    PROJECT_NAME,
    RECORD_IDENTIFIER,
    SAMPLE_NAME,
    SCHEMA_DESC_KEY,
    SCHEMA_ITEMS_KEY,
    SCHEMA_PROP_KEY,
    SCHEMA_TYPE_KEY,
    STATUS,
)
from pipestat.exceptions import PipestatError, SchemaError
from pipestat.parsed_schema import ParsedSchema

_LOGGER = logging.getLogger(__name__)


NULL_MAPPING_VALUE = {}
SCHEMA_PIPELINE_NAME_KEY = "pipeline_name"


# The columns associated with the file and image types
PATH_COL_SPEC = (Path, ...)
TITLE_COL_SPEC = (Optional[str], Field(default=None))
THUMBNAIL_COL_SPEC = (Optional[Path], Field(default=None))


def _custom_types_column_specifications():
    """Collection of the column specifications for the custom types"""
    return {
        "path": PATH_COL_SPEC,
        "title": TITLE_COL_SPEC,
        "thumbnail": THUMBNAIL_COL_SPEC,
    }


def get_base_model():
    class BaseModel(SQLModel):
        __table_args__ = {"extend_existing": True}

        model_config = ConfigDict(arbitrary_types_allowed=True)

    # return SQLModel
    return BaseModel


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


class ParsedSchemaDB(ParsedSchema):
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

    # def __init__(self, data: Union[Dict[str, Any], Path, str]) -> None:
    #     # initial validation and parse
    # def __init__(self):

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

    def _make_field_definitions(self, data: Dict[str, Any], require_type: bool):
        # TODO: default to string if no type key?
        # TODO: parse "required" ?
        defs = {}
        for name, subdata in data.items():
            result_indexed = False
            try:
                typename = subdata[SCHEMA_TYPE_KEY]
            except KeyError:
                if require_type:
                    _LOGGER.error(f"'{SCHEMA_TYPE_KEY}' is required for each schema element")
                    raise
                else:
                    data_type = str
            else:
                data_type = self._get_data_type(typename)
            if data_type == CLASSES_BY_TYPE["object"] or data_type == CLASSES_BY_TYPE["array"]:
                if "index" in subdata and subdata["index"] is True:
                    _LOGGER.warning(f"Cannot index JSONB Column, ignoring index: True for {name} ")
                defs[name] = (
                    data_type,
                    Field(
                        sa_column=Column(JSONB),
                        default=null(),
                    ),
                )
            else:
                if "index" in subdata:
                    if isinstance(subdata["index"], bool):
                        result_indexed = subdata["index"]
                defs[name] = (
                    data_type,
                    Field(
                        default=subdata.get("default"),
                        nullable=True,
                        index=result_indexed,
                    ),
                )
        return defs

    @staticmethod
    def _get_data_type(type_name):
        t = CLASSES_BY_TYPE[type_name]
        # return ARRAY if t == list else t
        return t

    @property
    def file_like_table_name(self):
        return self._table_name("files")

    def build_history_model(self, pipeline_type):
        """Creates model for history ORM
        :param str pipeline_type: project or sample-level pipeline
        :return model: (model, table_name)
        """
        if pipeline_type == "project":
            history_table_name = self.project_table_name + "_history"
            data = self.project_level_data
            main_table_id = self.project_table_name + ".id"

        elif pipeline_type == "sample":
            history_table_name = self.sample_table_name + "_history"
            data = self.sample_level_data
            main_table_id = self.sample_table_name + ".id"

        else:
            raise PipestatError(
                f"Building model requires pipeline type. Provided type: '{pipeline_type}' "
            )

        if not self.sample_level_data and not self.project_level_data:
            return None

        field_defs = self._make_field_definitions(data, require_type=True)
        field_defs = self._add_status_field(field_defs)
        field_defs = self._add_record_identifier_field(field_defs)
        field_defs = self._add_id_field(field_defs)
        field_defs = self._add_pipeline_name_field(field_defs)
        field_defs = self._add_created_time_field(field_defs)
        field_defs = self._add_modified_time_field(field_defs)

        field_defs["source_record_id"] = (
            int,
            Field(
                default=None,
                foreign_key=main_table_id,
            ),
        )

        history_model = _create_model(history_table_name, **field_defs)

        return history_model, history_table_name

    def build_model(self, pipeline_type):
        if pipeline_type == "project":
            data = self.project_level_data
            # if using the same output schema and thus, pipeline name for samples and project
            # we must ensure there are distinct table names in the same database.
            table_name = self.project_table_name

        elif pipeline_type == "sample":
            data = self.sample_level_data
            table_name = self.sample_table_name
        else:
            raise PipestatError(
                f"Building model requires pipeline type. Provided type: '{pipeline_type}' "
            )

        if not self.sample_level_data and not self.project_level_data:
            return None

        field_defs = self._make_field_definitions(data, require_type=True)
        field_defs = self._add_status_field(field_defs)
        field_defs = self._add_record_identifier_field(field_defs)
        field_defs = self._add_id_field(field_defs)
        # field_defs = self._add_project_name_field(field_defs)
        field_defs = self._add_pipeline_name_field(field_defs)
        field_defs = self._add_created_time_field(field_defs)
        field_defs = self._add_modified_time_field(field_defs)
        return _create_model(table_name, **field_defs)

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

    @staticmethod
    def _add_project_name_field(field_defs: Dict[str, Any]) -> Dict[str, Any]:
        if PROJECT_NAME in field_defs:
            raise SchemaError(
                f"'{PROJECT_NAME}' is reserved as identifier and can't be part of schema."
            )

        field_defs[PROJECT_NAME] = (str, Field(default=None, nullable=True))

        return field_defs

    @staticmethod
    def _add_pipeline_name_field(field_defs: Dict[str, Any]) -> Dict[str, Any]:
        if PIPELINE_NAME in field_defs:
            raise SchemaError(
                f"'{PIPELINE_NAME}' is reserved as identifier and can't be part of schema."
            )

        field_defs[PIPELINE_NAME] = (str, Field(default=None, nullable=True))

        return field_defs

    @staticmethod
    def _add_id_field(field_defs: Dict[str, Any]) -> Dict[str, Any]:
        if ID_KEY in field_defs:
            raise SchemaError(
                f"'{ID_KEY}' is reserved for primary key and can't be part of schema."
            )
        field_defs[ID_KEY] = (
            Optional[int],
            Field(default=None, primary_key=True),
        )
        return field_defs

    @staticmethod
    def _add_record_identifier_field(field_defs: Dict[str, Any]) -> Dict[str, Any]:
        if RECORD_IDENTIFIER in field_defs:
            raise SchemaError(
                f"'{RECORD_IDENTIFIER}' is reserved as identifier and can't be part of schema."
            )

        field_defs[RECORD_IDENTIFIER] = (str, Field(default=None, nullable=True))

        return field_defs

    @staticmethod
    def _add_sample_name_field(field_defs: Dict[str, Any]) -> Dict[str, Any]:
        if SAMPLE_NAME in field_defs:
            raise SchemaError(
                f"'{SAMPLE_NAME}' is reserved as identifier and can't be part of schema."
            )

        field_defs[SAMPLE_NAME] = (str, Field(default=None, nullable=True))

        return field_defs

    @staticmethod
    def _add_status_field(field_defs: Dict[str, Any]) -> Dict[str, Any]:
        if STATUS in field_defs:
            raise SchemaError(
                f"'{STATUS}' is reserved for status reporting and can't be part of schema."
            )

        field_defs[STATUS] = (str, Field(default=None, nullable=True))

        return field_defs

    @staticmethod
    def _add_created_time_field(field_defs: Dict[str, Any]) -> Dict[str, Any]:
        if CREATED_TIME in field_defs:
            raise SchemaError(
                f"'{CREATED_TIME}' is reserved for time reporting and can't be part of schema."
            )

        field_defs[CREATED_TIME] = (datetime.datetime, Field(default=None, nullable=True))

        return field_defs

    @staticmethod
    def _add_modified_time_field(field_defs: Dict[str, Any]) -> Dict[str, Any]:
        if MODIFIED_TIME in field_defs:
            raise SchemaError(
                f"'{MODIFIED_TIME}' is reserved for time reporting and can't be part of schema."
            )

        field_defs[MODIFIED_TIME] = (datetime.datetime, Field(default=None, nullable=True))

        return field_defs

    def _table_name(self, suffix: str) -> str:
        return f"{self.pipeline_name}__{suffix}"


def _create_model(table_name: str, **kwargs):
    return create_model(
        table_name,
        __base__=get_base_model(),
        __cls_kwargs__={"table": True},
        **kwargs,
    )


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
