import datetime
import os
from abc import ABC
from collections.abc import MutableMapping
from copy import deepcopy
from logging import getLogger
from typing import Any, Dict, Iterator, List, Optional, Union

from jsonschema import validate
from ubiquerg import mkabs
from yacman import FutureYAMLConfigManager as YAMLConfigManager
from yacman import load_yaml
from yacman.yacman_future import select_config

from pipestat.backends.file_backend.filebackend import FileBackend

from .const import (
    CFG_DATABASE_KEY,
    CFG_SCHEMA,
    CONFIG_KEY,
    CREATED_TIME,
    DATA_KEY,
    DB_ONLY_KEY,
    DB_URL,
    DEFAULT_PIPELINE_NAME,
    ENV_VARS,
    FILE_KEY,
    MODIFIED_TIME,
    MULTI_PIPELINE,
    OUTPUT_DIR,
    PIPELINE_NAME,
    PIPELINE_TYPE,
    PKG_NAME,
    PROJECT_NAME,
    RECORD_IDENTIFIER,
    RESULT_FORMATTER,
    SAMPLE_NAME_ID_KEY,
    SCHEMA_KEY,
    SCHEMA_PATH,
    STATUS_FILE_DIR,
    STATUS_SCHEMA,
    STATUS_SCHEMA_KEY,
    STATUS_SCHEMA_SOURCE_KEY,
)
from .exceptions import (
    ColumnNotFoundError,
    InvalidTimeFormatError,
    NoBackendSpecifiedError,
    PipestatDatabaseError,
    PipestatDependencyError,
    RecordNotFoundError,
    SchemaNotFoundError,
    PipestatSummarizeError,
)
from .helpers import default_formatter, make_subdirectories, validate_type, zip_report
from .reports import HTMLReportBuilder, _create_stats_objs_summaries

try:
    from pipestat.backends.db_backend.db_parsed_schema import ParsedSchemaDB as ParsedSchema
except ImportError:
    from .parsed_schema import ParsedSchema

try:
    from pipestat.backends.db_backend.db_helpers import construct_db_url
    from pipestat.backends.db_backend.dbbackend import DBBackend
except ImportError:
    # We let this pass, but if the user attempts to create DBBackend, check_dependencies raises exception.
    pass

try:
    from pipestat.backends.pephub_backend.pephubbackend import PEPHUBBACKEND
except ImportError:
    # Let this pass, if phc dependencies cannot be imported, raise exception
    pass


_LOGGER = getLogger(PKG_NAME)


def check_dependencies(dependency_list: list = None, msg: str = None):
    """Decorator to check that the dependency list has successfully been imported."""

    def wrapper(func):
        def inner(*args, **kwargs):
            dependencies_satisfied = True
            if dependency_list is not None:
                for i in dependency_list:
                    if i not in globals():
                        _LOGGER.warning(msg=f"Missing dependency: {i}")
                        dependencies_satisfied = False
            if dependencies_satisfied is False:
                raise PipestatDependencyError(msg=msg)
            return func(*args, **kwargs)

        return inner

    return wrapper


def require_backend(func):
    """Decorator to ensure a backend exists for functions that require one"""

    def inner(self, *args, **kwargs):
        if not self.backend:
            raise NoBackendSpecifiedError
        return func(self, *args, **kwargs)

    return inner


class PipestatManager(MutableMapping):
    """
    Pipestat standardizes reporting of pipeline results and
    pipeline status management. It formalizes a way for pipeline developers
    and downstream tools developers to communicate -- results produced by a
    pipeline can easily and reliably become an input for downstream analyses.
    A PipestatManager object exposes an API for interacting with the results and
    pipeline status and can be backed by either a YAML-formatted file
    or a database.
    """

    def __init__(
        self,
        project_name: Optional[str] = None,
        record_identifier: Optional[str] = None,
        schema_path: Optional[str] = None,
        results_file_path: Optional[str] = None,
        database_only: Optional[bool] = True,
        config_file: Optional[str] = None,
        config_dict: Optional[dict] = None,
        flag_file_dir: Optional[str] = None,
        show_db_logs: bool = False,
        pipeline_type: Optional[str] = None,
        pipeline_name: Optional[str] = None,
        result_formatter: staticmethod = default_formatter,
        multi_pipelines: bool = False,
        output_dir: Optional[str] = None,
        pephub_path: Optional[str] = None,
    ):
        """
        Initialize the PipestatManager object

        :param str record_identifier: record identifier to report for. This
            creates a weak bound to the record, which can be overridden in
            this object method calls
        :param str schema_path: path to the output schema that formalizes
            the results structure
        :param str results_file_path: YAML file to report into, if file is
            used as the object back-end
        :param bool database_only: whether the reported data should not be
            stored in the memory, but only in the database
        :param str config_file: path to the configuration file
        :param dict config_dict:  a mapping with the config file content
        :param str flag_file_dir: path to directory containing flag files
        :param bool show_db_logs: Defaults to False, toggles showing database logs
        :param str pipeline_type: "sample" or "project"
        :param str pipeline_name: name of the current pipeline, defaults to
        :param str result_formatter: function for formatting result
        :param bool multi_pipelines: allows for running multiple pipelines for one file backend
        :param str output_dir: target directory for report generation via summarize and table generation via table.
        """

        super(PipestatManager, self).__init__()

        # Initialize the cfg dict as an attribute that holds all configuration data
        self.cfg = {}

        # Load and validate database configuration
        # If results_file_path exists, backend is a file else backend is database.

        self.cfg["config_path"] = select_config(config_file, ENV_VARS["config"])

        if config_dict is not None:
            self.cfg[CONFIG_KEY] = YAMLConfigManager.from_obj(entries=config_dict)
        elif self.cfg["config_path"] is not None:
            self.cfg[CONFIG_KEY] = YAMLConfigManager.from_yaml_file(
                filepath=self.cfg["config_path"]
            )
        else:
            self.cfg[CONFIG_KEY] = YAMLConfigManager()

        cfg_schema = load_yaml(CFG_SCHEMA)
        validate(self.cfg[CONFIG_KEY].exp, cfg_schema)

        self.cfg["pephub_path"] = self.cfg[CONFIG_KEY].priority_get(
            "pephub_path", override=pephub_path
        )

        self.cfg[SCHEMA_PATH] = self.cfg[CONFIG_KEY].priority_get(
            "schema_path", env_var=ENV_VARS["schema"], override=schema_path
        )
        self.process_schema(schema_path)

        self.cfg[RECORD_IDENTIFIER] = self.cfg[CONFIG_KEY].priority_get(
            "record_identifier", env_var=ENV_VARS["record_identifier"], override=record_identifier
        )

        # TODO this is a work around for Looper ~ https://github.com/pepkit/looper/issues/492, sharing pipeline names
        # In the future, we should get pipeline name only from output schema.
        if pipeline_name:
            self.cfg[PIPELINE_NAME] = pipeline_name
        elif self.cfg[CONFIG_KEY].get(PIPELINE_NAME):
            self.cfg[PIPELINE_NAME] = self.cfg[CONFIG_KEY].get(PIPELINE_NAME)
        elif self.cfg[SCHEMA_KEY] and self.cfg[SCHEMA_KEY].pipeline_name:
            self.cfg[PIPELINE_NAME] = self.cfg[SCHEMA_KEY].pipeline_name
        else:
            self.cfg[PIPELINE_NAME] = DEFAULT_PIPELINE_NAME

        self.cfg[PROJECT_NAME] = self.cfg[CONFIG_KEY].priority_get(
            "project_name", env_var=ENV_VARS["project_name"], override=project_name
        )

        self.cfg[SAMPLE_NAME_ID_KEY] = self.cfg[CONFIG_KEY].priority_get(
            "record_identifier",
            env_var=ENV_VARS["sample_name"],
            override=record_identifier,
        )

        self.cfg[DB_ONLY_KEY] = database_only

        self.cfg[PIPELINE_TYPE] = self.cfg[CONFIG_KEY].priority_get(
            "pipeline_type", default="sample", override=pipeline_type
        )

        self.cfg[FILE_KEY] = mkabs(
            self.resolve_results_file_path(
                self.cfg[CONFIG_KEY].priority_get(
                    "results_file_path",
                    env_var=ENV_VARS["results_file"],
                    override=results_file_path,
                )
            ),
            self.cfg["config_path"],
        )

        if "{record_identifier}" in str(self.cfg[FILE_KEY]):
            # In the special case where the user wants to use {record_identifier} in file path
            pass
        else:
            make_subdirectories(self.cfg[FILE_KEY])

        self.cfg[RESULT_FORMATTER] = result_formatter

        self.cfg[MULTI_PIPELINE] = multi_pipelines

        self.cfg["multi_result_files"] = None

        self.cfg[OUTPUT_DIR] = self.cfg[CONFIG_KEY].priority_get("output_dir", override=output_dir)

        if self.cfg[FILE_KEY]:
            self.initialize_filebackend(record_identifier, results_file_path, flag_file_dir)

        elif self.cfg["pephub_path"]:
            self.initialize_pephubbackend(record_identifier, self.cfg["pephub_path"])
        else:
            self.initialize_dbbackend(record_identifier, show_db_logs)

    def __str__(self):
        """
        Generate string representation of the object

        :return str: string representation of the object
        """
        res = f"{self.__class__.__name__} ({self.cfg[PIPELINE_NAME]})"
        res += "\nBackend: {}".format(
            f"File\n - results: {self.cfg[FILE_KEY]}\n - status: {self.cfg[STATUS_FILE_DIR]}"
            if self.cfg[FILE_KEY]
            else f"Database (dialect: {self.backend.db_engine_key})"
        )
        if self.cfg[FILE_KEY]:
            res += f"\nMultiple Pipelines Allowed: {self.cfg[MULTI_PIPELINE]}"
        else:
            res += f"\nProject Name: {self.cfg[PROJECT_NAME]}"
            res += f"\nDatabase URL: {self.cfg[DB_URL]}"
            res += f"\nConfig File: {self.config_path}"

        res += f"\nPipeline name: {self.cfg[PIPELINE_NAME]}"
        res += f"\nPipeline type: {self.cfg[PIPELINE_TYPE]}"
        if self.cfg[SCHEMA_PATH] is not None:
            res += "\nProject Level Data:"
            for k, v in self.cfg[SCHEMA_KEY].project_level_data.items():
                res += f"\n {k} : {v}"
            res += "\nSample Level Data:"
            for k, v in self.cfg[SCHEMA_KEY].sample_level_data.items():
                res += f"\n {k} : {v}"
        res += f"\nStatus Schema key: {self.cfg[STATUS_SCHEMA_KEY]}"
        res += f"\nResults formatter: {str(self.cfg[RESULT_FORMATTER].__name__)}"
        res += f"\nResults schema source: {self.cfg[SCHEMA_PATH]}"
        res += f"\nStatus schema source: {self.cfg[STATUS_SCHEMA_SOURCE_KEY]}"
        res += f"\nRecords count: {self.record_count}"
        if self.cfg[SCHEMA_PATH] is not None:
            high_res = self.highlighted_results
        else:
            high_res = None
        if high_res:
            res += f"\nHighlighted results: {', '.join(high_res)}"
        return res

    def __getitem__(self, key):
        # This is a wrapper for the retrieve function:
        result = self.retrieve_one(record_identifier=key)
        return result

    def __setitem__(self, key, value):
        # This is a wrapper for the report function:
        result = self.report(record_identifier=key, values=value)
        return result

    def __delitem__(self, key):
        # This is a wrapper for the remove function; it removes the entire record:
        result = self.remove(record_identifier=key)
        return result

    def __iter__(
        self,
        limit: Optional[int] = 1000,
        cursor: Optional[int] = None,
    ) -> Iterator:
        """
        This is a wrapper around select_records that creates an iterator of records

        :param int limit: maximum number of results to retrieve per page
        :param int cursor: cursor position to begin retrieving records
        :return: Iterator
        """
        if self.file:
            # File backend does not support cursor-based paging
            return iter(self.select_records(limit=limit)["records"])
        else:
            return iter(self.select_records(limit=limit, cursor=cursor)["records"])

    def __len__(self):
        return len(self.cfg)

    def resolve_results_file_path(self, results_file_path):
        """Replace {record_identifier} in results_file_path if it exists.
        :param str results_file_path: YAML file to report into, if file is
        used as the object back-end
        """
        # Save for later when assessing if there may be multiple result files
        if results_file_path:
            assert isinstance(results_file_path, str), TypeError("Path is expected to be a str")
            if self.record_identifier:
                try:
                    self.cfg["unresolved_result_path"] = results_file_path
                    results_file_path = results_file_path.format(
                        record_identifier=self.record_identifier
                    )
                    return results_file_path
                except AttributeError:
                    self.cfg["unresolved_result_path"] = results_file_path
                    return results_file_path
            else:
                self.cfg["unresolved_result_path"] = results_file_path
                return results_file_path
        return results_file_path

    def initialize_filebackend(
        self,
        record_identifier: str = None,
        results_file_path: str = None,
        flag_file_dir: str = None,
    ):
        """
        Initializes the file backend
        :param str record_identifier: the record identifier
        :param str results_file_path: the path to the results file used for the backend
        :param str flag_file_dir: the path to the flag file directory
        """

        # Check if there will be multiple results_file_paths
        _LOGGER.debug(f"Determined file as backend: {results_file_path}")

        if self.cfg[DB_ONLY_KEY]:
            _LOGGER.debug(
                "Running in database only mode does not make sense with a YAML file as a backend. "
                "Changing back to using memory."
            )
            self.cfg[DB_ONLY_KEY] = False

        flag_file_dir = self.cfg[CONFIG_KEY].priority_get(
            "flag_file_dir",
            override=flag_file_dir,
            default=os.path.dirname(self.cfg[FILE_KEY]),
        )
        self.cfg[STATUS_FILE_DIR] = mkabs(flag_file_dir, self.config_path or self.cfg[FILE_KEY])
        make_subdirectories(self.cfg[STATUS_FILE_DIR])

        self.backend = FileBackend(
            self.cfg[FILE_KEY],
            record_identifier,
            self.cfg[PIPELINE_NAME],
            self.cfg[PIPELINE_TYPE],
            self.cfg[SCHEMA_KEY],
            self.cfg[STATUS_SCHEMA_KEY],
            self.cfg[STATUS_FILE_DIR],
            self.cfg[RESULT_FORMATTER],
            self.cfg[MULTI_PIPELINE],
        )

        return

    def initialize_pephubbackend(self, record_identifier: str = None, pephub_path: str = None):
        """
        Initializes the pephub backend
        :param str record_identifier: the record identifier
        :param str pephub_path: the path to the pephub registry
        """
        self.backend = PEPHUBBACKEND(
            record_identifier,
            pephub_path,
            self.cfg[PIPELINE_NAME],
            self.cfg[PIPELINE_TYPE],
            self.cfg[SCHEMA_KEY],
            self.cfg[STATUS_SCHEMA_KEY],
            self.cfg[RESULT_FORMATTER],
        )

    @check_dependencies(
        dependency_list=["DBBackend"],
        msg="Missing required dependencies for this usage, e.g. try pip install pipestat['dbbackend']",
    )
    def initialize_dbbackend(self, record_identifier: str = None, show_db_logs: bool = False):
        """
        Initializes the database backend
        :param str record_identifier: the record identifier
        :param bool show_db_logs: boolean to show_db_logs
        """
        _LOGGER.debug("Determined database as backend")
        if self.cfg[SCHEMA_KEY] is None:
            raise SchemaNotFoundError("Output schema must be supplied for DB backends.")
        if CFG_DATABASE_KEY not in self.cfg[CONFIG_KEY]:
            raise NoBackendSpecifiedError()
        try:
            dbconf = self.cfg[CONFIG_KEY].exp[
                CFG_DATABASE_KEY
            ]  # the .exp expands the paths before url construction
            if "sqlite_url" in dbconf:
                sqlite_url = f"sqlite:///{dbconf['sqlite_url']}"
                self.cfg[DB_URL] = sqlite_url
            else:
                self.cfg[DB_URL] = construct_db_url(dbconf)
        except KeyError:
            raise PipestatDatabaseError(f"No database section ('{CFG_DATABASE_KEY}') in config")
        self._show_db_logs = show_db_logs

        self.backend = DBBackend(
            record_identifier,
            self.cfg[PIPELINE_NAME],
            show_db_logs,
            self.cfg[PIPELINE_TYPE],
            self.cfg[SCHEMA_KEY],
            self.cfg[STATUS_SCHEMA_KEY],
            self.cfg[DB_URL],
            self.cfg[STATUS_SCHEMA_SOURCE_KEY],
            self.cfg[RESULT_FORMATTER],
        )

    @require_backend
    def clear_status(
        self,
        record_identifier: str = None,
        flag_names: List[str] = None,
    ) -> List[Union[str, None]]:
        """
        Remove status flags

        :param str record_identifier: name of the sample_level record to remove flags for
        :param Iterable[str] flag_names: Names of flags to remove, optional; if
            unspecified, all schema-defined flag names will be used.
        :return List[Union[str,None]]: Collection of names of flags removed
        """

        r_id = record_identifier or self.record_identifier
        return self.backend.clear_status(record_identifier=r_id, flag_names=flag_names)

    @require_backend
    def count_records(self) -> int:
        """
        Count records
        :return int: number of records
        """
        return self.backend.count_records()

    @require_backend
    def get_status(
        self,
        record_identifier: str = None,
    ) -> Optional[str]:
        """
        Get the current pipeline status
        :param str record_identifier: name of the sample_level record
        :return str: status identifier, e.g. 'running'

        """

        r_id = record_identifier or self.record_identifier
        return self.backend.get_status(record_identifier=r_id)

    @require_backend
    def list_recent_results(
        self,
        limit: Optional[int] = 1000,
        start: Optional[str] = None,
        end: Optional[str] = None,
        time_column: Optional[str] = "modified",
    ) -> dict:
        """
        :param int  limit: limit number of results returned
        :param str start: most recent result to filter on, defaults to now, e.g. 2023-10-16 13:03:04.680400
        :param str end: oldest result to filter on, e.g. 1970-10-16 13:03:04.680400
        :param str time_column: created or modified column/attribute to filter on
        :return dict results: a dict containing start, end, num of records, and list of retrieved records
        """

        if self.cfg["pephub_path"]:
            _LOGGER.warning(f"List recent results not supported for PEPHub backend")
            return {}
        date_format = "%Y-%m-%d %H:%M:%S"
        if start is None:
            start = datetime.datetime.now()
        else:
            try:
                start = datetime.datetime.strptime(start, date_format)
            except ValueError:
                raise InvalidTimeFormatError(msg=f"Incorrect time format, requires:{date_format}")

        if end is None:
            end = datetime.datetime.strptime("1900-01-01 00:00:00", date_format)
        else:
            try:
                end = datetime.datetime.strptime(end, date_format)
            except ValueError:
                raise InvalidTimeFormatError(msg=f"Incorrect time format, requires: {date_format}")

        if time_column == "created":
            col_name = CREATED_TIME
        else:
            col_name = MODIFIED_TIME

        results = self.select_records(
            limit=limit,
            filter_conditions=[
                {
                    "key": col_name,
                    "operator": "lt",
                    "value": start,
                },
                {
                    "key": col_name,
                    "operator": "gt",
                    "value": end,
                },
            ],
        )
        return results

    def process_schema(self, schema_path):
        # Load pipestat schema in two parts: 1) main and 2) status
        self._schema_path = self.cfg[CONFIG_KEY].priority_get(
            "schema_path", env_var=ENV_VARS["schema"], override=schema_path
        )

        if self._schema_path is None:
            _LOGGER.warning("No pipestat output schema was supplied to PipestatManager.")
            self.cfg[SCHEMA_KEY] = None
            self.cfg[STATUS_SCHEMA_KEY] = None
            self.cfg[STATUS_SCHEMA_SOURCE_KEY] = None
            # return None
            # raise SchemaNotFoundError("PipestatManager creation failed; no schema")
        else:
            # Main schema
            schema_to_read = mkabs(self._schema_path, self.cfg["config_path"])
            self._schema_path = schema_to_read
            parsed_schema = ParsedSchema(schema_to_read)
            self.cfg[SCHEMA_KEY] = parsed_schema

            # Status schema
            self.cfg[STATUS_SCHEMA_KEY] = parsed_schema.status_data
            if not self.cfg[STATUS_SCHEMA_KEY]:
                self.cfg[STATUS_SCHEMA_SOURCE_KEY] = STATUS_SCHEMA
                self.cfg[STATUS_SCHEMA_KEY] = load_yaml(filepath=STATUS_SCHEMA)
            else:
                self.cfg[STATUS_SCHEMA_SOURCE_KEY] = schema_to_read

    @require_backend
    def remove(
        self,
        record_identifier: str = None,
        result_identifier: str = None,
    ) -> bool:
        """
        Remove a result.

        If no result ID specified or last result is removed, the entire record
        will be removed.

        :param str record_identifier: name of the sample_level record
        :param str result_identifier: name of the result to be removed or None
             if the record should be removed.
        :return bool: whether the result has been removed
        """

        r_id = record_identifier or self.cfg[RECORD_IDENTIFIER]
        return self.backend.remove(
            record_identifier=r_id,
            result_identifier=result_identifier,
        )

    @require_backend
    def remove_record(
        self,
        record_identifier: Optional[str] = None,
        rm_record: Optional[bool] = False,
    ) -> bool:
        return self.backend.remove_record(
            record_identifier=record_identifier,
            rm_record=rm_record,
        )

    @require_backend
    def report(
        self,
        values: Dict[str, Any],
        record_identifier: Optional[str] = None,
        force_overwrite: bool = True,
        result_formatter: Optional[staticmethod] = None,
        strict_type: bool = True,
        history_enabled: bool = True,
    ) -> Union[List[str], bool]:
        """
        Report a result.

        :param Dict[str, any] values: dictionary of result-value pairs
        :param str record_identifier: unique identifier of the record, value
            in 'record_identifier' column to look for to determine if the record
            already exists
        :param bool force_overwrite: whether to overwrite the existing record
        :param str result_formatter: function for formatting result
        :param bool strict_type: whether the type of the reported values should
            remain as is. Pipestat would attempt to convert to the
            schema-defined one otherwise
        :param bool history_enabled: Should history of reported results be enabled?
        :return str reported_results: return list of formatted string
        """

        result_formatter = result_formatter or self.cfg[RESULT_FORMATTER]
        values = deepcopy(values)
        r_id = record_identifier or self.cfg[RECORD_IDENTIFIER]
        if r_id is None:
            raise NotImplementedError("You must supply a record identifier to report results")

        result_identifiers = list(values.keys())
        if self.cfg[SCHEMA_KEY] is not None:
            for r in result_identifiers:
                # First confirm this property is defined in the schema
                if r not in self.result_schemas:
                    raise ColumnNotFoundError(
                        f"Can't report a result for attribute '{r}'; it is not defined in the output schema."
                    )

                validate_type(
                    value=values[r],
                    schema=self.result_schemas[r],
                    strict_type=strict_type,
                    record_identifier=record_identifier,
                )

        reported_results = self.backend.report(
            values=values,
            record_identifier=r_id,
            force_overwrite=force_overwrite,
            result_formatter=result_formatter,
            history_enabled=history_enabled,
        )

        return reported_results

    @require_backend
    def select_distinct(
        self,
        columns: Optional[Union[str, List[str]]] = None,
    ) -> List[Any]:
        """
        Retrieves unique results for a list of attributes.

        :param List[str] columns: columns to include in the result
        :return list[any] result: this is a list of distinct results
        """
        if not isinstance(columns, list) and not isinstance(columns, str):
            raise ValueError(
                "Columns must be a list of strings or string, e.g. ['record_identifier', 'number_of_things']"
            )

        result = self.backend.select_distinct(columns=columns)
        return result

    @require_backend
    def select_records(
        self,
        columns: Optional[List[str]] = None,
        filter_conditions: Optional[List[Dict[str, Any]]] = None,
        limit: Optional[int] = 1000,
        cursor: Optional[int] = None,
        bool_operator: Optional[str] = "AND",
    ) -> Dict[str, Any]:
        """
        :param list[str] columns: columns to include in the result
        :param list[dict]  filter_conditions: e.g. [{"key": ["id"], "operator": "eq", "value": 1)], operator list:
            - eq for ==
            - lt for <
            - ge for >=
            - in for in_
            - like for like
        :param int limit: maximum number of results to retrieve per page
        :param int cursor: cursor position to begin retrieving records
        :param bool bool_operator: Perform filtering with AND or OR Logic.
        :return dict records_dict = {
            "total_size": int,
            "page_size": int,
            "next_page_token": int,
            "records": List[Dict[{key, Any}]],
        }
        """

        return self.backend.select_records(
            columns=columns,
            filter_conditions=filter_conditions,
            limit=limit,
            cursor=cursor,
            bool_operator=bool_operator,
        )

    @require_backend
    def retrieve_one(
        self,
        record_identifier: str = None,
        result_identifier: Optional[Union[str, List[str]]] = None,
    ) -> Union[Any, Dict[str, Any]]:
        """
        Retrieve a single record
        :param str record_identifier: single record_identifier
        :param str result_identifier: single record_identifier or list of result identifiers
        :return: Dict[str, any]: a mapping with filtered results reported for the record
        """
        record_identifier = record_identifier or self.record_identifier

        filter_conditions = [
            {
                "key": "record_identifier",
                "operator": "eq",
                "value": record_identifier,
            },
        ]
        if result_identifier:
            if isinstance(result_identifier, str):
                columns = [result_identifier]
            elif isinstance(result_identifier, list):
                columns = result_identifier
            else:
                raise ValueError("Result identifier must be a str or list[str]")
            result = self.select_records(filter_conditions=filter_conditions, columns=columns)[
                "records"
            ]
            if len(result) > 0:
                if len(columns) > 1:
                    try:
                        return result[0]
                    except IndexError:
                        raise RecordNotFoundError(
                            f"Results '{columns}' for '{record_identifier}' not found"
                        )
                try:
                    return result[0][columns[0]]
                except IndexError:
                    raise RecordNotFoundError(
                        f"Results '{columns}' for '{record_identifier}' not found"
                    )
            else:
                raise RecordNotFoundError(
                    f"Results '{columns}' for '{record_identifier}' not found"
                )
        else:
            try:
                result = self.select_records(filter_conditions=filter_conditions)["records"]
                if len(result) > 0:
                    try:
                        return result[0]
                    except IndexError:
                        raise RecordNotFoundError(f"Record '{record_identifier}' not found")
                else:
                    raise RecordNotFoundError(f"Record '{record_identifier}' not found")
            except IndexError:
                raise RecordNotFoundError(f"Record '{record_identifier}' not found")

    def retrieve_history(
        self,
        record_identifier: str = None,
        result_identifier: Optional[Union[str, List[str]]] = None,
    ) -> Union[Any, Dict[str, Any]]:
        """
        Retrieve a single record's history
        :param str record_identifier: single record_identifier
        :param str result_identifier: single result_identifier or list of result identifiers
        :return: Dict[str, any]: a mapping with filtered historical results
        """

        record_identifier = record_identifier or self.record_identifier

        if self.file:
            result = self.backend.select_records(
                filter_conditions=[
                    {
                        "key": "record_identifier",
                        "operator": "eq",
                        "value": record_identifier,
                    }
                ],
                meta_data_bool=True,
            )["records"][0]

            if "meta" in result and "history" in result["meta"]:
                history = {}
                if isinstance(result_identifier, str) and result_identifier in result:
                    history.update(
                        {result_identifier: result["meta"]["history"][result_identifier]}
                    )
                elif isinstance(result_identifier, list):
                    for r in result_identifier:
                        if r in result["meta"]["history"]:
                            history.update({r: result["meta"]["history"][r]})
                else:
                    history = result["meta"]["history"]
            else:
                _LOGGER.warning(f"No history available for Record: {record_identifier}")
                return {}

        elif self.cfg["pephub_path"]:
            _LOGGER.warning(f"Retrieving history not supported for PEPHub backend")
            return None
        else:
            if result_identifier:
                history = self.backend.retrieve_history_db(record_identifier, result_identifier)[
                    "history"
                ]
            else:
                history = self.backend.retrieve_history_db(record_identifier)["history"]

            # DB backend returns some extra_keys that we can remove before returning them to the user.
            extra_keys_to_delete = [
                "id",
                "pipestat_created_time",
                "source_record_id",
                "record_identifier",
            ]
            history = {
                key: value for key, value in history.items() if key not in extra_keys_to_delete
            }

        return history

    def retrieve_many(
        self,
        record_identifiers: List[str],
        result_identifier: Optional[str] = None,
    ) -> Union[Any, Dict[str, Any]]:
        """
        :param record_identifiers: list of record identifiers
        :param str result_identifier: single record_identifier
        :return: Dict[str, any]: a mapping with filtered
            results reported for the record
        """

        filter = {
            "key": "record_identifier",
            "operator": "in",
            "value": record_identifiers,
        }
        if result_identifier:
            result = self.select_records(filter_conditions=[filter], columns=[result_identifier])
        else:
            result = self.select_records(filter_conditions=[filter])

        if len(result["records"]) == 0:
            RecordNotFoundError(f"Records, '{record_identifiers}',  not found")
        else:
            return result

    @require_backend
    def set_status(
        self,
        status_identifier: str,
        record_identifier: str = None,
    ) -> None:
        """
        Set pipeline run status.

        The status identifier needs to match one of identifiers specified in
        the status schema. A basic, ready to use, status schema is shipped with
        this package.

        :param str status_identifier: status to set, one of statuses defined
            in the status schema
        :param str record_identifier: sample_level record identifier to set the
            pipeline status for
        """
        r_id = record_identifier or self.record_identifier
        self.backend.set_status(status_identifier, r_id)

    @require_backend
    def link(self, link_dir) -> Union[str, None]:
        """
        This function creates a link structure such that results are organized by type.
        :param str link_dir: path to desired symlink output directory
        :return str | None linked_results_path: path to symlink directory or None
        """

        self.check_multi_results()
        linked_results_path = self.backend.link(link_dir=link_dir)

        return linked_results_path

    @require_backend
    def summarize(
        self,
        looper_samples: Optional[list] = None,
        amendment: Optional[str] = None,
        portable: Optional[bool] = False,
        output_dir: Optional[str] = None,
    ) -> Union[str, None]:
        """
        Builds a browsable html report for reported results.
        :param Iterable[str] looper_samples: list of looper Samples from PEP
        :param Iterable[str] amendment: name indicating amendment to use, optional
        :param bool portable: moves figures and report files to directory for easy sharing
        :param str output_dir: overrides output_dir set during pipestatManager creation.
        :return str: report_path

        """

        if output_dir:
            self.cfg[OUTPUT_DIR] = output_dir

        if self.cfg["pephub_path"]:
            if OUTPUT_DIR not in self.cfg:
                _LOGGER.warning(f"Output directory is required for pipestat summarize.")
                return None

        self.check_multi_results()

        # Before proceeding check if there are any results at the specified backend
        try:
            current_results = self.select_records()
            if len(current_results["records"]) < 1:
                raise PipestatSummarizeError(f"No results found at specified backend")
        except Exception as e:
            raise PipestatSummarizeError(f"PipestatSummarizeError due to exception: {e}")

        html_report_builder = HTMLReportBuilder(prj=self, portable=portable)
        report_path = html_report_builder(
            pipeline_name=self.cfg[PIPELINE_NAME],
            amendment=amendment,
            looper_samples=looper_samples,
        )

        if portable is True:
            report_path = zip_report(report_dir_name=os.path.dirname(report_path))

        return report_path

    def check_multi_results(self):
        # Check to see if the user used a path with "{record-identifier}"
        if self.file:
            if "{record_identifier}" in self.cfg["unresolved_result_path"]:
                # assume there are multiple result files in sub-directories
                self.cfg["multi_result_files"] = True
                results_directory = self.cfg["unresolved_result_path"].split(
                    "{record_identifier}"
                )[0]
                results_directory = mkabs(results_directory, self.cfg["config_path"])
                make_subdirectories(results_directory)
                self.backend.aggregate_multi_results(results_directory)

    @require_backend
    def table(
        self,
        output_dir: Optional[str] = None,
    ) -> List[str]:
        """
        Generates stats (.tsv) and object (.yaml) files.
        :param str output_dir: overrides output_dir set during pipestatManager creation.
        :return list[str] table_path_list: list containing output file paths of stats and objects

        """
        if output_dir:
            self.cfg[OUTPUT_DIR] = output_dir

        self.check_multi_results()
        pipeline_name = self.cfg[PIPELINE_NAME]
        table_path_list = _create_stats_objs_summaries(self, pipeline_name)

        return table_path_list

    def _get_attr(self, attr: str) -> Any:
        """
        Safely get the name of the selected attribute of this object

        :param str attr: attr to select
        :return:
        """
        return self.get(attr)

    @property
    def config_path(self) -> str:
        """
        Config path. None if the config was not provided or if provided
        as a mapping of the config contents

        :return str: path to the provided config
        """
        return self.cfg.get("config_path", None)

    @property
    def data(self) -> YAMLConfigManager:
        """
        Data object

        :return yacman.YAMLConfigManager: the object that stores the reported data
        """
        return self.backend._data

    @property
    def db_url(self) -> str:
        """
        Database URL, generated based on config credentials

        :return str: database URL
        :raise PipestatDatabaseError: if the object is not backed by a database
        """
        return self.cfg[DB_URL]

    @property
    def file(self) -> str:
        """
        File path that the object is reporting the results into

        :return str: file path that the object is reporting the results into
        """
        return self.cfg[FILE_KEY]

    @property
    def highlighted_results(self) -> List[str]:
        """
        Highlighted results

        :return List[str]: a collection of highlighted results
        """
        return [k for k, v in self.result_schemas.items() if v.get("highlight") is True]

    @property
    def output_dir(self) -> str:
        """
        Output directory for report and stats generation

        :return str: path to output_dir
        """
        return self.cfg[OUTPUT_DIR]

    @property
    def pipeline_name(self) -> str:
        """
        Pipeline name

        :return str: Pipeline name
        """
        return self.cfg[PIPELINE_NAME]

    @property
    def project_name(self) -> str:
        """
        Project name the object writes the results to

        :return str: project name the object writes the results to
        """
        return self.cfg[PROJECT_NAME]

    @property
    def pipeline_type(self) -> str:
        """
        Pipeline type: "sample" or "project"

        :return str: pipeline type
        """
        return self.cfg[PIPELINE_TYPE]

    @property
    def record_identifier(self) -> str:
        """
        Pipeline type: "sample" or "project"

        :return str: pipeline type
        """
        return self.cfg[RECORD_IDENTIFIER]

    @property
    def record_count(self) -> int:
        """
        Number of records reported

        :return int: number of records reported
        """
        return self.count_records()

    @property
    def result_schemas(self) -> Dict[str, Any]:
        """
        Result schema mappings

        :return dict: schemas that formalize the structure of each result
            in a canonical jsonschema way
        """
        return {
            **self.cfg[SCHEMA_KEY].project_level_data,
            **self.cfg[SCHEMA_KEY].sample_level_data,
        }

    @property
    def schema(self) -> ParsedSchema:
        """
        Schema mapping

        :return ParsedSchema: schema object that formalizes the results structure
        """
        return self.cfg["_schema"]

    @property
    def schema_path(self) -> str:
        """
        Schema path

        :return str: path to the provided schema
        """
        return self.cfg[SCHEMA_PATH]

    @property
    def status_schema(self) -> Dict:
        """
        Status schema mapping

        :return dict: schema that formalizes the pipeline status structure
        """
        return self.cfg[STATUS_SCHEMA_KEY]

    @property
    def status_schema_source(self) -> Dict:
        """
        Status schema source

        :return dict: source of the schema that formalizes
            the pipeline status structure
        """
        return self.cfg[STATUS_SCHEMA_SOURCE_KEY]


class SamplePipestatManager(PipestatManager):
    def __init__(self, **kwargs):
        PipestatManager.__init__(self, pipeline_type="sample", **kwargs)
        _LOGGER.warning("Initialize PipestatMgrSample")


class ProjectPipestatManager(PipestatManager):
    def __init__(self, **kwargs):
        PipestatManager.__init__(self, pipeline_type="project", **kwargs)
        _LOGGER.warning("Initialize PipestatMgrProject")


class PipestatBoss(ABC):
    """
    PipestatBoss simply holds Sample or Project Managers that are child classes of PipestatManager.
        :param list[str] pipeline_list: list that holds pipeline types, e.g. ['sample','project']
        :param str record_identifier: record identifier to report for. This
            creates a weak bound to the record, which can be overridden in
            this object method calls
        :param str schema_path: path to the output schema that formalizes
            the results structure
        :param str results_file_path: YAML file to report into, if file is
            used as the object back-end
        :param bool database_only: whether the reported data should not be
            stored in the memory, but only in the database
        :param str | dict config: path to the configuration file or a mapping
            with the config file content
        :param str flag_file_dir: path to directory containing flag files
        :param bool show_db_logs: Defaults to False, toggles showing database logs
        :param str pipeline_type: "sample" or "project"
        :param str result_formatter: function for formatting result
        :param bool multi_pipelines: allows for running multiple pipelines for one file backend
        :param str output_dir: target directory for report generation via summarize and table generation via table.
    """

    def __init__(self, pipeline_list: Optional[list] = None, **kwargs):
        _LOGGER.warning("Initialize PipestatBoss")
        if len(pipeline_list) > 3:
            _LOGGER.warning(
                "PipestatBoss currently only supports one 'sample' and one 'project' pipeline. Ignoring extra types."
            )
        for i in pipeline_list:
            if i == "sample":
                self.samplemanager = SamplePipestatManager(**kwargs)
            elif i == "project":
                self.projectmanager = ProjectPipestatManager(**kwargs)
            else:
                _LOGGER.warning(f"This pipeline type is not supported. Pipeline supplied: {i}")

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)
