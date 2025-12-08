import datetime
import operator
import os.path
from copy import deepcopy
from functools import reduce
from glob import glob
from itertools import chain
from logging import getLogger
from typing import Any, Callable, Dict, List, Literal, Optional, Tuple, Union

from ubiquerg import create_lock, remove_lock
from yacman import FutureYAMLConfigManager as YAMLConfigManager
from yacman import read_lock, write_lock

from ...backends.abstract import PipestatBackend
from ...const import CREATED_TIME, DATE_FORMAT, HISTORY_KEY, META_KEY, MODIFIED_TIME, PKG_NAME
from ...exceptions import PipestatError, UnrecognizedStatusError
from ...helpers import get_all_result_files

_LOGGER = getLogger(PKG_NAME)


class FileBackend(PipestatBackend):
    def __init__(
        self,
        results_file_path: str,
        record_identifier: Optional[str] = None,
        pipeline_name: Optional[str] = None,
        pipeline_type: Optional[str] = None,
        parsed_schema: Optional[str] = None,
        status_schema: Optional[str] = None,
        status_file_dir: Optional[str] = None,
        result_formatter: Optional[staticmethod] = None,
        multi_pipelines: Optional[bool] = None,
    ):
        """
        Class representing a File backend
        :param str results_file_path: YAML file to report into, if file is
            used as the object back-end
        :param str record_identifier: record identifier to report for. This
            creates a weak bound to the record, which can be overridden in
            this object method calls
        :param str pipeline_name: name of pipeline associated with result
        :param str pipeline_type: "sample" or "project"
        :param str parsed_schema: results output schema.
        :param str status_schema: schema containing pipeline statuses e.g. 'running'
        :param str status_file_dir: directory for placing status flags
        :param str result_formatter: function for formatting result
        :param bool multi_pipelines: allows for running multiple pipelines for one file backend

        """
        super().__init__(pipeline_type)
        _LOGGER.debug("Initialize FileBackend")

        self.results_file_path = results_file_path
        self.pipeline_name = pipeline_name
        self.pipeline_type = pipeline_type
        self.record_identifier = record_identifier
        self.parsed_schema = parsed_schema
        self.status_schema = status_schema
        self.status_file_dir = status_file_dir
        self.result_formatter = result_formatter
        self.multi_pipelines = multi_pipelines

        self.determine_results_file()

    def determine_results_file(self) -> None:
        """Initialize or load results_file from given path
        :param str results_file_path: YAML file to report into, if file is
        used as the object back-end
        """

        if "{record_identifier}" in self.results_file_path:
            # In the special case where the user wants to use {record_identifier} in file path
            pass
        else:
            if not os.path.exists(self.results_file_path):
                _LOGGER.debug(
                    f"Results file doesn't yet exist. Initializing: {self.results_file_path}"
                )
                self._init_results_file()
            else:
                _LOGGER.debug(f"Loading results file: {self.results_file_path}")
                self._load_results_file()

    def check_record_exists(
        self,
        record_identifier: str,
    ) -> bool:
        """
        Check if the specified record exists in self._data

        :param str record_identifier: record to check for
        :return bool: whether the record exists in the table
        """

        return (
            self.pipeline_name in self._data
            and record_identifier in self._data[self.pipeline_name][self.pipeline_type]
        )

    def clear_status(
        self, record_identifier: str = None, flag_names: List[str] = None
    ) -> List[Union[str, None]]:
        """
        Remove status flags

        :param str record_identifier: name of the record to remove flags for
        :param Iterable[str] flag_names: Names of flags to remove, optional; if
            unspecified, all schema-defined flag names will be used.
        :return List[str]: Collection of names of flags removed
        """

        flag_names = flag_names or list(self.status_schema.keys())
        if isinstance(flag_names, str):
            flag_names = [flag_names]
        removed = []
        for f in flag_names:
            path_flag_file = self.get_status_flag_path(
                status_identifier=f, record_identifier=record_identifier
            )
            try:
                os.remove(path_flag_file)
            except:
                pass
            else:
                _LOGGER.info(f"Removed existing flag: {path_flag_file}")
                removed.append(f)
        return removed

    def count_records(self):
        """
        Count records
        :return int: number of records
        """

        return len(self._data[self.pipeline_name])

    def get_flag_file(self, record_identifier: str = None) -> Union[str, List[str], None]:
        """
        Get path to the status flag file for the specified record

        :param str record_identifier: unique record identifier
        :return str | list[str] | None: path to the status flag file
        """

        r_id = record_identifier
        regex = os.path.join(self.status_file_dir, f"{self.pipeline_name}_{r_id}_*.flag")
        file_list = glob(regex)
        if len(file_list) > 1:
            _LOGGER.warning("Multiple flag files found")
            return file_list
        elif len(file_list) == 1:
            return file_list[0]
        else:
            _LOGGER.debug("No flag files found")
            return None
        pass

    def get_status(self, record_identifier: str) -> Optional[str]:
        """
        Get the current pipeline status

        :param str record_identifier: record identifier to set the
            pipeline status for
        :return str: status identifier, e.g. 'running'
        """
        r_id = record_identifier or self.record_identifier
        flag_file = self.get_flag_file(record_identifier=record_identifier)
        if flag_file is not None:
            assert isinstance(flag_file, str), TypeError(
                "Flag file path is expected to be a str, were multiple flags found?"
            )
            with open(flag_file, "r") as f:
                status = f.read()
            return status
        _LOGGER.debug(
            f"Could not determine status for '{r_id}' record. "
            f"No flags found in: {self.status_file_dir}"
        )
        return None

    def get_status_flag_path(self, status_identifier: str, record_identifier=None) -> str:
        """
        Get the path to the status file flag

        :param str status_identifier: one of the defined status IDs in schema
        :param str record_identifier: unique record ID
        :return str: absolute path to the flag file or None if object is
            backed by a DB
        """

        r_id = record_identifier
        return os.path.join(
            self.status_file_dir,
            f"{self.pipeline_name}_{r_id}_{status_identifier}.flag",
        )

    def list_results(
        self,
        restrict_to: Optional[List[str]] = None,
        record_identifier: Optional[str] = None,
    ) -> List[str]:
        """
        Lists all, or a selected set of, reported results

        :param List[str] restrict_to: selected subset of names of results to list
        :param str record_identifier: unique identifier of the record
        :return List[str]: names of results which exist
        """
        record_identifier = record_identifier or self.record_identifier

        try:
            results = list(
                self._data[self.pipeline_name][self.pipeline_type][record_identifier].keys()
            )
        except KeyError:
            return []
        if restrict_to:
            return [r for r in restrict_to if r in results]
        return results

    def remove(
        self,
        record_identifier: Optional[str] = None,
        result_identifier: Optional[str] = None,
    ) -> bool:
        """
        Remove a result.

        If no result ID specified or last result is removed, the entire record
        will be removed.

        :param str record_identifier: unique identifier of the record
        :param str result_identifier: name of the result to be removed or None
             if the record should be removed.
        :return bool: whether the result has been removed
        """

        record_identifier = record_identifier or self.record_identifier

        rm_record = True if result_identifier is None else False

        if not self.check_record_exists(
            record_identifier=record_identifier,
        ):
            _LOGGER.error(f"Record '{record_identifier}' not found")
            return False

        if result_identifier and not self.check_result_exists(
            result_identifier, record_identifier
        ):
            _LOGGER.error(f"'{result_identifier}' has not been reported for '{record_identifier}'")
            return False

        if rm_record:
            # NOTE: THIS CURRENTLY REMOVES ALL HISTORY OF THE RECORD AS WELL
            self.remove_record(
                record_identifier=record_identifier,
                rm_record=rm_record,
            )
        else:
            self._modify_history(
                data=self._data[self.pipeline_name][self.pipeline_type][record_identifier],
                res_id=result_identifier,
                time=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                value="",
            )
            del self._data[self.pipeline_name][self.pipeline_type][record_identifier][
                result_identifier
            ]
            _LOGGER.info(
                f"Removed result '{result_identifier}' for record "
                f"'{record_identifier}' from '{self.pipeline_name}' namespace"
            )
            if not self._data[self.pipeline_name][self.pipeline_type][record_identifier]:
                _LOGGER.info(
                    f"Last result removed for '{record_identifier}'. " f"Removing the record"
                )
                rm_record = True
                self.remove_record(
                    record_identifier=record_identifier,
                    rm_record=rm_record,
                )
            # Check if the last remaining attributes are the timestamps
            remaining_attributes = list(
                self._data[self.pipeline_name][self.pipeline_type][record_identifier].keys()
            )
            if len(remaining_attributes) == 1 and META_KEY in remaining_attributes:
                _LOGGER.info(f"Last result removed for '{record_identifier}'.")

            with write_lock(self._data) as locked_data:
                locked_data.write()
        return True

    def remove_record(
        self,
        record_identifier: Optional[str] = None,
        rm_record: Optional[bool] = False,
    ) -> bool:
        """
        Remove a record, requires rm_record to be True

        :param str record_identifier: unique identifier of the record
        :param bool rm_record: bool for removing record.
        :return bool: whether the result has been removed
        """
        if rm_record:
            try:
                _LOGGER.info(f"Removing '{record_identifier}' record")
                del self._data[self.pipeline_name][self.pipeline_type][record_identifier]
                with write_lock(self._data) as data_locked:
                    data_locked.write()
                return True
            except:
                _LOGGER.warning(
                    f" Unable to remove record, aborting Removing '{record_identifier}' record"
                )
                return False
        else:
            _LOGGER.info(
                f" rm_record flag is set to False, aborting Removing '{record_identifier}' record"
            )

    def report(
        self,
        values: Dict[str, Any],
        record_identifier: Optional[str] = None,
        force_overwrite: bool = True,
        result_formatter: Optional[staticmethod] = None,
        history_enabled: bool = True,
    ) -> Union[List[str], bool]:
        """
        Update the value of a result in a current namespace.

        This method overwrites any existing data and creates the required
         hierarchical mapping structure if needed.

        :param Dict[str, Any] values: dict of results identifiers and values
            to be reported
        :param str record_identifier: unique identifier of the record
        :param bool force_overwrite: Toggles force overwriting results, defaults to False
        :param str result_formatter: function for formatting result
        :return bool | list[str] results_formatted: return list of formatted string
        """

        # record_identifier = record_identifier or self.record_identifier
        record_identifier = record_identifier

        result_formatter = result_formatter or self.result_formatter
        results_formatted = []
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        result_identifiers = list(values.keys())
        if self.parsed_schema is not None:
            self.assert_results_defined(
                results=result_identifiers, pipeline_type=self.pipeline_type
            )

        existing = self.list_results(
            record_identifier=record_identifier,
            restrict_to=result_identifiers,
        )

        if existing:
            existing_str = ", ".join(existing)
            _LOGGER.warning(f"These results exist for '{record_identifier}': {existing_str}")
            if not force_overwrite:
                return False
            _LOGGER.info(f"Overwriting existing results: {existing_str}")
            self._data[self.pipeline_name][self.pipeline_type][record_identifier][META_KEY].update(
                {MODIFIED_TIME: current_time}
            )
        if not existing:
            if record_identifier in self._data[self.pipeline_name][self.pipeline_type].keys():
                self._data[self.pipeline_name][self.pipeline_type][record_identifier].setdefault(
                    META_KEY, {}
                )
                self._data[self.pipeline_name][self.pipeline_type][record_identifier][
                    META_KEY
                ].update({MODIFIED_TIME: current_time})

            else:
                self._data[self.pipeline_name][self.pipeline_type].setdefault(
                    record_identifier, {}
                )
                self._data[self.pipeline_name][self.pipeline_type][record_identifier].setdefault(
                    META_KEY, {}
                )
                self._data[self.pipeline_name][self.pipeline_type][record_identifier][
                    META_KEY
                ].update({MODIFIED_TIME: current_time})
                self._data[self.pipeline_name][self.pipeline_type][record_identifier][
                    META_KEY
                ].update({CREATED_TIME: current_time})

        for res_id, val in values.items():
            if history_enabled:
                if existing:
                    self._modify_history(
                        data=self._data[self.pipeline_name][self.pipeline_type][record_identifier][
                            META_KEY
                        ],
                        res_id=res_id,
                        time=current_time,
                        value=self._data[self.pipeline_name][self.pipeline_type][
                            record_identifier
                        ][res_id],
                    )
            self._data[self.pipeline_name][self.pipeline_type][record_identifier][res_id] = val
            results_formatted.append(
                result_formatter(
                    pipeline_name=self.pipeline_name,
                    record_identifier=record_identifier,
                    res_id=res_id,
                    value=val,
                )
            )

        with write_lock(self._data) as locked_data:
            locked_data.write()

        return results_formatted

    def select_distinct(self, columns: Union[str, List[str]]) -> List[Tuple]:
        """
        Retrieve distinct results given a list of columns

        :param str | List[str] columns: column columns to include in the result
        :return List[Tuple]: returns distinct values.
        """
        if isinstance(columns, str):
            columns = [columns]

        record_data = self._data.data[self.pipeline_name][self.pipeline_type]

        final_values_list = []
        for record_id in record_data.keys():
            record_value_list = []
            for column in columns:
                if column in list(record_data[record_id].keys()):
                    record_value_list.append(record_data[record_id][column])
                else:
                    if column != "record_identifier":
                        _LOGGER.warning(msg=f"Column {column} not reported, skipping....")

            if "record_identifier" in columns:
                # record_identifier is not a column in the file backend but the user may want it
                record_value_list.append(record_id)
            final_values_list.append(tuple(record_value_list))

        # make sure that the order is correct
        unique_lists = list(set(final_values_list))
        return unique_lists

    def select_records(
        self,
        columns: Optional[List[str]] = None,
        filter_conditions: Optional[List[Dict[str, Any]]] = None,
        limit: Optional[int] = 1000,
        cursor: Optional[int] = None,
        bool_operator: Optional[str] = "AND",
        meta_data_bool: Optional[bool] = False,
    ) -> Dict[str, Any]:
        """
        Select records from the FileBackend

        :param list[str] columns: columns to include in the result
        :param list[dict]  filter_conditions: e.g. [{"key": ["id"], "operator": "eq", "value": 1)], operator list:
            - eq for ==
            - lt for <
            - ge for >=
            - in for in_
        :param int limit: maximum number of results to retrieve per page
        :param int cursor: cursor position to begin retrieving records
        :param bool bool_operator: Perform filtering with AND or OR Logic.
        :param bool meta_data: Should this return associated meta data with records?
        :return dict records_dict = {
            "total_size": int,
            "page_size": int,
            "next_page_token": int,
            "records": List[Dict[{key, Any}]],
        }
        """
        if cursor:
            _LOGGER.warning("Cursor not supported for FileBackend, ignoring cursor")

        def get_operator(op: Literal["eq", "lt", "ge", "gt", "in"]) -> Any:
            """
            Get python operator for a given string

            :param str op: desired operator, "eq", "lt"
            :return: operator function
            """

            if op == "eq":
                return operator.__eq__
            if op == "lt":
                return operator.__lt__
            if op == "ge":
                return operator.__ge__
            if op == "gt":
                return operator.__gt__
            if op == "in":
                return operator.contains
            raise ValueError(f"Invalid filter operator: {op}")

        def get_nested_column(result_value: dict, key_list: list, retrieved_operator: Callable):
            """
            Recursive function that evaluates a nested list of keys vs a value dict

            :param dict result_value: nested dictionary
            :param key_list: list of keys, e.g. keys that may be in the nested dictionary e.g. ['id', 'name']
            :param retrieved_operator: operator function (ge, gt...)
            :return bool:
            """
            if len(key_list) == 1:
                if result_value.get(key_list[0], None):
                    if retrieved_operator(key_list, result_value.get(key_list[0])):
                        return True
                return False
            else:
                return get_nested_column(
                    result_value[key_list[0]], key_list[1:], retrieved_operator
                )

        records_list = []

        data = deepcopy(self._data.data[self.pipeline_name][self.pipeline_type])

        if columns:
            if not isinstance(columns, list):
                raise ValueError(
                    "Columns must be a list of strings, e.g. ['record_identifier', 'number_of_things']"
                )

        total_count = len(data.keys())

        filtered_records_list = []
        if filter_conditions:
            for filter_condition in filter_conditions:
                if list(filter_condition.keys()) != ["key", "operator", "value"]:
                    raise ValueError(
                        "Filter conditions must be a dictionary with keys 'key', 'operator', and 'value'"
                    )

                retrieved_operator = get_operator(filter_condition["operator"])
                retrieved_results = []

                # Check each sample's dict
                for record_identifier in list(data.keys())[0:limit]:
                    if filter_condition["key"] != "record_identifier":
                        for key, value in data[record_identifier].items():
                            result = False
                            if isinstance(value, dict) and key != "meta":
                                if key == filter_condition["key"][0]:
                                    result = get_nested_column(
                                        value,
                                        filter_condition["key"][1:],
                                        retrieved_operator,
                                    )
                            else:
                                if key == "meta":
                                    # Filter datetime objects
                                    if (
                                        filter_condition["key"] == CREATED_TIME
                                        or filter_condition["key"] == MODIFIED_TIME
                                    ):
                                        try:
                                            time_stamp = datetime.datetime.strptime(
                                                data[record_identifier][META_KEY][
                                                    filter_condition["key"]
                                                ],
                                                DATE_FORMAT,
                                            )
                                            result = retrieved_operator(
                                                time_stamp, filter_condition["value"]
                                            )
                                        except TypeError:
                                            result = False
                                elif filter_condition["key"] == key:
                                    result = retrieved_operator(value, filter_condition["value"])

                            if result:
                                retrieved_results.append(record_identifier)
                    else:
                        # If user wants record_identifier
                        if isinstance(filter_condition["value"], list):
                            for v in filter_condition["value"]:
                                if (
                                    record_identifier == v
                                    and record_identifier not in retrieved_results
                                ):
                                    retrieved_results.append(record_identifier)
                        elif record_identifier == filter_condition["value"]:
                            retrieved_results.append(record_identifier)

                if retrieved_results:
                    filtered_records_list.append(retrieved_results)
        else:
            # Assume user wants all the records if no filter was given.
            filtered_records_list = [list(data.keys())[0:limit]]

        # There is now a list of dicts for each filtered condition.
        # Depending on Union or Intersection we want to pare down the list.
        shared_keys = []

        if bool_operator.lower() == "and" and filtered_records_list:
            shared_keys = list(reduce(lambda i, j: i & j, (set(x) for x in filtered_records_list)))

        if bool_operator.lower() == "or" and filtered_records_list:
            shared_keys = list(set(chain(*filtered_records_list)))

        if shared_keys:
            for record_identifier in sorted(shared_keys):
                record = {}
                if columns:  # Did the user specify a list of columns as well?
                    for key, value in list(data[record_identifier].items()):
                        if key in columns:
                            record.update({key: value})
                    if "record_identifier" in columns:
                        record.update({"record_identifier": record_identifier})
                else:
                    record = data[record_identifier]
                if record != {}:
                    record.update({"record_identifier": record_identifier})
                    records_list.append(record)
                    if "meta" in record and not meta_data_bool:
                        del record["meta"]

        records_dict = {
            "total_size": total_count,
            "page_size": limit,
            "next_page_token": 0,
            "records": records_list,
        }

        return records_dict

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
        :param str record_identifier: record identifier to set the
            pipeline status for
        """
        r_id = record_identifier or self.record_identifier
        if self.status_schema is not None:
            known_status_identifiers = self.status_schema.keys()
            if status_identifier not in known_status_identifiers:
                raise UnrecognizedStatusError(
                    f"'{status_identifier}' is not a defined status identifier. "
                    f"These are allowed: {known_status_identifiers}"
                )
        prev_status = self.get_status(r_id)

        if prev_status is not None:
            prev_flag_path = self.get_status_flag_path(prev_status, record_identifier)
            os.remove(prev_flag_path)
        flag_path = self.get_status_flag_path(status_identifier, record_identifier)
        create_lock(flag_path)
        with open(flag_path, "w") as f:
            f.write(status_identifier)
        remove_lock(flag_path)

        if prev_status:
            _LOGGER.debug(f"Changed status from '{prev_status}' to '{status_identifier}'")

    def _htmlreportbuilder(self):
        """
        build html report based on all reported results
        """

        # build new folder for the report
        self.reports_dir = os.path.join(self.results_file_path, "reports")
        _LOGGER.debug(f"Reports dir: {self.reports_dir}")

    def _init_results_file(self) -> None:
        """
        Initialize YAML results file if it does not exist.
        Read the data stored in the existing file into the memory otherwise.

        :return bool: whether the file has been created
        """
        _LOGGER.info(f"Initializing results file '{self.results_file_path}'")

        # Must ensure sub-directories exist if they do not
        # TODO should this actually be handled by yacman?
        try:
            os.makedirs(os.path.dirname(self.results_file_path))
        except FileExistsError:
            pass

        self._data = YAMLConfigManager.from_yaml_file(
            self.results_file_path,
            create_file=True,
        )
        self._data.update_from_obj({self.pipeline_name: {}})
        self._data.setdefault(self.pipeline_name, {})
        self._data[self.pipeline_name].setdefault("project", {})
        self._data[self.pipeline_name].setdefault("sample", {})
        with write_lock(self._data) as data_locked:
            data_locked.write()

    def aggregate_multi_results(self, results_directory) -> None:
        """
        Collects single results files and aggregates them into a new aggregate_results.yaml file
        :param str results_directory: directory containing subdirectories containing results.yaml files
        :return: None
        """
        all_result_files = get_all_result_files(results_directory)
        if len(all_result_files) == 0:
            _LOGGER.warning(
                "Attempting to aggregate multiple results files but no result files found. Ensure they are in subdirectories, e.g. 'record1/record1_results.yaml'"
            )
        aggregate_results_file_path = os.path.join(results_directory, "aggregate_results.yaml")

        # THIS WILL OVERWRITE self.results_file_path and self._data on the current psm!
        self.results_file_path = aggregate_results_file_path
        self._init_results_file()

        for file in all_result_files:
            try:
                temp_data = YAMLConfigManager.from_yaml_file(file)
            except ValueError:
                temp_data = YAMLConfigManager()
            if self.pipeline_name in temp_data:
                if "project" in temp_data[self.pipeline_name]:
                    self._data[self.pipeline_name]["project"].update(
                        temp_data[self.pipeline_name]["project"]
                    )
                if "sample" in temp_data[self.pipeline_name]:
                    self._data[self.pipeline_name]["sample"].update(
                        temp_data[self.pipeline_name]["sample"]
                    )

        with write_lock(self._data) as data_locked:
            data_locked.write()

    def _load_results_file(self) -> None:
        _LOGGER.debug(f"Reading data from '{self.results_file_path}'")
        data = YAMLConfigManager.from_yaml_file(self.results_file_path)
        if not bool(data):
            self._data = data
            self._data.setdefault(self.pipeline_name, {})
            self._data[self.pipeline_name].setdefault("project", {})
            self._data[self.pipeline_name].setdefault("sample", {})
            with write_lock(self._data) as data_locked:
                data_locked.write()
        namespaces_reported = [k for k in data.keys() if not k.startswith("_")]
        num_namespaces = len(namespaces_reported)
        if num_namespaces == 0:
            self._data = data
        elif num_namespaces > 0:
            if self.pipeline_name in namespaces_reported:
                self._data = data
            elif self.pipeline_name not in namespaces_reported and self.multi_pipelines is True:
                self._data = data
                self._data.setdefault(self.pipeline_name, {})
                self._data[self.pipeline_name].setdefault("project", {})
                self._data[self.pipeline_name].setdefault("sample", {})
                _LOGGER.warning("MULTI PIPELINES FOR SINGLE RESULTS FILE")
            else:
                raise PipestatError(
                    f"Trying to report result for namespace '{self.pipeline_name}' at '{self.results_file_path}', but "
                    f"{num_namespaces} other namespaces are already in the file: [{', '.join(namespaces_reported)}]. "
                    f"Pipestat will not report multiple namespaces to one file unless `multi_pipelines` is True."
                )

    def _modify_history(self, data, res_id, time, value):
        """Modify File backend with each change

        data is the loaded yaml results file in dict format
        type = "report", "deletion"
        """
        if "history" not in data:
            data.setdefault(HISTORY_KEY, {})
        if res_id not in data[HISTORY_KEY]:
            data[HISTORY_KEY].setdefault(res_id, {})

        data[HISTORY_KEY][res_id].update({time: value})
