import copy
import datetime
from contextlib import contextmanager
from logging import getLogger
from typing import Any, Dict, List, NoReturn, Optional, Tuple, Union

from sqlmodel import Session, SQLModel, create_engine
from sqlmodel import select as sql_select

from pipestat.backends.abstract import PipestatBackend
from pipestat.backends.db_backend.db_helpers import selection_filter

from ...const import CREATED_TIME, MODIFIED_TIME, PKG_NAME, RECORD_IDENTIFIER, STATUS
from ...exceptions import (
    ColumnNotFoundError,
    PipestatDatabaseError,
    RecordNotFoundError,
    SchemaError,
    SchemaNotFoundError,
    UnrecognizedStatusError,
)

_LOGGER = getLogger(PKG_NAME)


class DBBackend(PipestatBackend):
    def __init__(
        self,
        record_identifier: Optional[str] = None,
        pipeline_name: Optional[str] = None,
        show_db_logs: bool = False,
        pipeline_type: Optional[str] = False,
        parsed_schema: Optional[str] = None,
        status_schema: Optional[str] = None,
        db_url: Optional[str] = None,
        status_schema_source: Optional[dict] = None,
        result_formatter: Optional[staticmethod] = None,
    ):
        """
        Class representing a Database backend
        :param str record_identifier: record identifier to report for. This
            creates a weak bound to the record, which can be overridden in
            this object method calls
        :param str pipeline_name: name of pipeline associated with result
        :param str show_db_logs: Defaults to False, toggles showing database logs
        :param str pipeline_type: "sample" or "project"
        :param str parsed_schema: results output schema. Used to construct DB columns.
        :param str status_schema: schema containing pipeline statuses e.g. 'running'
        :param str db_url: url used for connection to Postgres DB
        :param dict status_schema_source: filepath of status schema
        :param str result_formatter: function for formatting result
        """

        super().__init__(pipeline_type)
        _LOGGER.debug(f"Initializing DBBackend for pipeline '{pipeline_name}'")
        self.pipeline_name = pipeline_name
        self.pipeline_type = pipeline_type or "sample"
        self.record_identifier = record_identifier
        self.parsed_schema = parsed_schema
        self.status_schema = status_schema
        self.db_url = db_url
        self.show_db_logs = show_db_logs
        self.status_schema_source = status_schema_source
        self.result_formatter = result_formatter

        self.orms = self._create_orms(pipeline_type=pipeline_type)
        self.history_table = self._create_history_orms(pipeline_type=pipeline_type)

        self.table_name = list(self.orms.keys())[0]
        SQLModel.metadata.create_all(self._engine)

    def check_record_exists(
        self,
        record_identifier: str,
    ) -> bool:
        """
        Check if the specified record exists in the table

        :param str record_identifier: record to check for
        :return bool: whether the record exists in the table
        """

        query_hit = self.select_records(
            filter_conditions=[
                {
                    "key": "record_identifier",
                    "operator": "eq",
                    "value": record_identifier,
                }
            ]
        )

        return bool(query_hit["records"])

    def count_records(self):
        """
        Count rows in a selected table
        :return int: number of records
        """

        mod = self.get_model(table_name=self.table_name)
        with self.session as s:
            stmt = sql_select(mod)
            records = s.exec(stmt).all()
            return len(records)

    def get_model(self, table_name: str):
        """
        Get model based on table_name
        :param str table_name: pipelinename__sample or pipelinename__project
        :return mod: model/orm associated with the table name
        """
        if self.orms is None:
            raise PipestatDatabaseError("Object relational mapper classes not defined.")

        mod = self.orms.get(table_name)

        if mod is None:
            raise PipestatDatabaseError(
                f"No object relational mapper class defined for table '{table_name}'. "
                f"{len(self.orms)} defined: {', '.join(self.orms.keys())}"
            )
        return mod

    def get_status(self, record_identifier: str) -> Optional[str]:
        """
        Get pipeline status

        :param str record_identifier: record identifier to set the
            pipeline status for
        :return str status
        """

        try:
            result = self.select_records(
                columns=[STATUS],
                filter_conditions=[
                    {
                        "key": "record_identifier",
                        "operator": "eq",
                        "value": record_identifier,
                    }
                ],
            )
        except RecordNotFoundError:
            return None
        try:
            status = result["records"][0]["status"]
        except IndexError:
            status = None

        return status

    def list_results(
        self,
        restrict_to: Optional[List[str]] = None,
        record_identifier: str = None,
    ) -> List[str]:
        """
        Check if the specified results exist in the table

        :param List[str] restrict_to: results identifiers to check for
        :param str record_identifier: record to check for
        :return List[str] existing: if no result identifier specified, return all results for the record
        :return List[str]: results identifiers that exist
        """

        rid = record_identifier
        record = self.select_records(
            filter_conditions=[
                {
                    "key": "record_identifier",
                    "operator": "eq",
                    "value": rid,
                }
            ]
        )
        try:
            record = record["records"][0]
        except IndexError:
            return []

        if restrict_to is None:
            return (
                [
                    key
                    for key in self.parsed_schema.results_data.keys()
                    if getattr(record, key, None) is not None
                ]
                if record
                else []
            )
        else:
            return [r for r in restrict_to if record.get(r, None) is not None] if record else []

    def remove(
        self,
        record_identifier: Optional[str] = None,
        result_identifier: Optional[str] = None,
    ) -> bool:
        """
        Remove a result.

        If no result ID specified, the entire record
        will be removed.

        :param str record_identifier: unique identifier of the record
        :param str result_identifier: name of the result to be removed or None
             if the record should be removed.
        :return bool: whether the result has been removed
        """

        # TODO removing last result identifier (apart from created, modified time should remove record)
        record_identifier = record_identifier or self.record_identifier
        rm_record = True if result_identifier is None else False

        if not self.check_record_exists(record_identifier=record_identifier):
            _LOGGER.error(f"Record '{record_identifier}' not found")
            return False

        if result_identifier and not self.check_result_exists(
            result_identifier=result_identifier,
            record_identifier=record_identifier,
        ):
            _LOGGER.error(f"'{result_identifier}' has not been reported for '{record_identifier}'")
            return False

        try:
            ORMClass = self.get_model(table_name=self.table_name)
            if self.check_record_exists(
                record_identifier=record_identifier,
            ):
                with self.session as s:
                    records = s.exec(
                        sql_select(ORMClass).where(
                            getattr(ORMClass, "record_identifier") == record_identifier
                        )
                    )
                    if rm_record is True:
                        self.remove_record(
                            record_identifier=record_identifier,
                            rm_record=rm_record,
                        )
                    else:
                        if not self.check_result_exists(
                            record_identifier=record_identifier,
                            result_identifier=result_identifier,
                        ):
                            raise RecordNotFoundError(
                                f"Result '{result_identifier}' not found for record "
                                f"'{record_identifier}'"
                            )
                        setattr(records.first(), result_identifier, None)
                    s.commit()
            else:
                raise RecordNotFoundError(f"Record '{record_identifier}' not found")
        except Exception as e:
            _LOGGER.error(f"Could not remove the result from the database. Exception: {e}")
            raise

        return True

    def remove_record(
        self,
        record_identifier: Optional[str] = None,
        rm_record: Optional[bool] = False,
    ) -> NoReturn:
        """
        Remove a record, requires rm_record to be True

        :param str record_identifier: unique identifier of the record
        :param bool rm_record: bool for removing record.
        :return bool: whether the result has been removed
        :raises RecordNotFoundError: if record not found
        """

        record_identifier = record_identifier or self.record_identifier
        if rm_record:
            try:
                ORMClass = self.get_model(table_name=self.table_name)
                ORMClass_History = self.history_table[list(self.history_table.keys())[0]]
                if self.check_record_exists(
                    record_identifier=record_identifier,
                ):
                    with self.session as s:
                        source_record_id = (
                            s.exec(
                                sql_select(ORMClass).where(
                                    getattr(ORMClass, RECORD_IDENTIFIER) == record_identifier
                                )
                            )
                            .first()
                            .id
                        )
                        linked_records = s.exec(
                            sql_select(ORMClass_History).where(
                                getattr(ORMClass_History, "source_record_id") == source_record_id
                            )
                        ).all()
                        for r in linked_records:
                            s.delete(r)
                        s.commit()
                    with self.session as s:
                        record = s.exec(
                            sql_select(ORMClass).where(
                                getattr(ORMClass, "record_identifier") == record_identifier
                            )
                        ).first()
                        s.delete(record)
                        s.commit()
                else:
                    raise RecordNotFoundError(f"Record '{record_identifier}' not found")
            except Exception as e:
                _LOGGER.error(f"Could not remove the result from the database. Exception: {e}")
                raise
        else:
            _LOGGER.info(f" rm_record flag False, aborting Removing '{record_identifier}' record")

    def report(
        self,
        values: Dict[str, Any],
        record_identifier: str,
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
        :param bool force_overwrite: force overwriting of results, defaults to False.
        :param str result_formatter: function for formatting result
        :return list[str] list results_formatted | bool: return list of formatted string
        """

        record_identifier = record_identifier or self.record_identifier

        result_formatter = result_formatter or self.result_formatter
        results_formatted = []
        result_identifiers = list(values.keys())
        if self.parsed_schema is None:
            raise SchemaNotFoundError("DB Backend report results requires schema")
        self.assert_results_defined(results=result_identifiers, pipeline_type=self.pipeline_type)

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

        try:
            ORMClass = self.get_model(table_name=self.table_name)
            ORMClass_History = self.history_table[list(self.history_table.keys())[0]]
            values.update({RECORD_IDENTIFIER: record_identifier})

            if not self.check_record_exists(
                record_identifier=record_identifier,
            ):
                current_time = datetime.datetime.now()
                values.update({CREATED_TIME: current_time})
                values.update({MODIFIED_TIME: current_time})
                new_record = ORMClass(**values)
                with self.session as s:
                    s.add(new_record)
                    s.commit()
            else:
                with self.session as s:
                    record_to_update = s.exec(
                        sql_select(ORMClass).where(
                            getattr(ORMClass, RECORD_IDENTIFIER) == record_identifier
                        )
                    ).first()
                    old_record_attributes = record_to_update.model_dump()
                    values.update({MODIFIED_TIME: datetime.datetime.now()})
                    for result_id, result_value in values.items():
                        setattr(record_to_update, result_id, result_value)
                    s.commit()
                if history_enabled:
                    if "id" in old_record_attributes:
                        del old_record_attributes["id"]
                    with self.session as s:
                        source_record = s.exec(
                            sql_select(ORMClass).where(
                                getattr(ORMClass, RECORD_IDENTIFIER) == record_identifier
                            )
                        ).first()
                        new_record_history = ORMClass_History(**old_record_attributes)
                        new_record_history.source_record_id = source_record.id
                        s.add(new_record_history)
                        s.commit()

            for res_id, val in values.items():
                results_formatted.append(
                    result_formatter(
                        pipeline_name=self.pipeline_name,
                        record_identifier=record_identifier,
                        res_id=res_id,
                        value=val,
                    )
                )
            return results_formatted
        except Exception as e:
            _LOGGER.error(f"Could not insert the result into the database. Exception: {e}")
            raise

    def select_records(
        self,
        columns: Optional[List[str]] = None,
        filter_conditions: Optional[List[Dict[str, Any]]] = None,
        limit: Optional[int] = 1000,
        cursor: Optional[int] = None,
        bool_operator: Optional[str] = "AND",
    ) -> Dict[str, Any]:
        """
        Perform a `SELECT` on the table

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

        ORM = self.get_model(table_name=self.table_name)

        with self.session as s:

            try:
                total_count = len(s.exec(sql_select(ORM)).all())
            except Exception as e:
                raise PipestatDatabaseError(
                    msg=f"Could not get total_count. Is the database empty? Original Error Message: {e}"
                )

            if columns is not None:
                columns = copy.deepcopy(columns)
                for i in ["id", "record_identifier"]:  # Must add id, need it for cursor
                    if i not in columns:
                        columns.insert(0, i)
                try:
                    statement = sql_select(*[getattr(ORM, column) for column in columns]).order_by(
                        ORM.id
                    )
                except AttributeError:
                    raise ColumnNotFoundError(
                        msg=f"One of the supplied columns does not exist in current table: {columns}"
                    )
            else:
                statement = sql_select(ORM).order_by(ORM.id)

            if cursor is not None:
                statement = statement.where(ORM.id > cursor)

            statement = selection_filter(
                ORM=ORM,
                statement=statement,
                filter_conditions=filter_conditions,
                bool_operator=bool_operator,
            )

            if isinstance(limit, int):
                statement = statement.limit(limit)

            results = s.exec(statement).all()

        if results != []:
            next_cursor = results[-1].id
        else:
            next_cursor = None

        end_results = []

        # SQL model returns either a SQLModelMetaCLass OR a sqlalchemy Row.
        # We must create a dictionary containing the record before returning
        if not columns:
            end_results = [r.model_dump() for r in results]

        else:
            for record in results:
                record_dict = dict(record._mapping)
                del record_dict["id"]
                end_results.append(record_dict)

        records_dict = {
            "total_size": total_count,
            "page_size": limit,
            "next_page_token": next_cursor,
            "records": end_results,
        }

        return records_dict

    def retrieve_history_db(
        self,
        record_identifier: str,
        result_identifier: Optional[Union[str, List[str]]] = None,
    ) -> Dict[str, Any]:
        """

        :param record_identifier: single record_identifier
        :param result_identifier: single or list of result identifiers
        :return: dict records_dict = {
            "history": List[Dict[{key, Any}]],
        }
        """

        record_identifier = record_identifier or self.record_identifier

        ORMClass = self.get_model(table_name=self.table_name)
        ORMClass_History = self.history_table[list(self.history_table.keys())[0]]

        if not result_identifier:
            columns = None
        else:
            if isinstance(result_identifier, str):
                columns = [result_identifier]
            elif isinstance(result_identifier, list):
                columns = copy.deepcopy(result_identifier)
            else:
                raise ValueError("Result identifier must be a str or list[str]")
            for i in ["id", MODIFIED_TIME]:
                if i not in columns:
                    columns.insert(0, i)

        if not self.check_record_exists(
            record_identifier=record_identifier,
        ):
            raise RecordNotFoundError(f"{record_identifier} does not exist.")
        else:
            with self.session as s:
                source_record_id = (
                    s.exec(
                        sql_select(ORMClass).where(
                            getattr(ORMClass, RECORD_IDENTIFIER) == record_identifier
                        )
                    )
                    .first()
                    .id
                )
                if columns is not None:
                    try:
                        statement = sql_select(
                            *[getattr(ORMClass_History, column) for column in columns]
                        ).order_by(ORMClass_History.id)
                    except AttributeError:
                        raise ColumnNotFoundError(
                            msg=f"One of the supplied columns does not exist in current table: {columns}"
                        )
                else:
                    statement = sql_select(ORMClass_History).order_by(ORMClass_History.id)

                statement = statement.where(
                    getattr(ORMClass_History, "source_record_id") == source_record_id
                )

                history_records = s.exec(statement).all()

        end_results = []

        # SQL model returns either a SQLModelMetaCLass OR a sqlalchemy Row.
        # We must create a dictionary containing the record before returning

        if not columns:
            end_results = [r.model_dump() for r in history_records]

        else:
            for record in history_records:
                record_dict = dict(record._mapping)
                end_results.append(record_dict)

        # This next step is to process the results such that they will match output similar to the filebackend

        collected_keys = []
        new_history_dict = {}
        for result in end_results:
            for key, value in result.items():
                if key == MODIFIED_TIME:
                    continue
                elif value:
                    if key not in new_history_dict:
                        collected_keys.append(key)
                        new_history_dict[key] = {result[MODIFIED_TIME]: value}
                    else:
                        new_history_dict[key].update({result[MODIFIED_TIME]: value})

        records_dict = {
            "history": new_history_dict,
        }

        return records_dict

    def select_distinct(
        self,
        columns: Union[str, List[str]],
    ) -> List[Tuple]:
        """
        Perform a `SELECT DISTINCT` on given table and column

        :param str | List[str] columns: columns to include in the result
        :return List[Tuple]: returns distinct values.
        """
        if isinstance(columns, str):
            columns = [columns]

        ORM = self.get_model(table_name=self.table_name)
        with self.session as s:
            list_columns = [getattr(ORM, column) for column in columns]
            result = s.exec(sql_select(*list_columns).distinct()).all()

        return result

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

        record_identifier = record_identifier or self.record_identifier
        known_status_identifiers = self.status_schema.keys()
        if status_identifier not in known_status_identifiers:
            raise UnrecognizedStatusError(
                f"'{status_identifier}' is not a defined status identifier. "
                f"These are allowed: {known_status_identifiers}"
            )
        prev_status = self.get_status(record_identifier)
        try:
            self.report(
                values={STATUS: status_identifier},
                record_identifier=record_identifier,
            )
        except Exception as e:
            _LOGGER.error(
                f"Could not insert into the status table ('{self.table_name}'). Exception: {e}"
            )
            raise
        if prev_status:
            _LOGGER.debug(f"Changed status from '{prev_status}' to '{status_identifier}'")

    def _create_orms(self, pipeline_type):
        """Create ORMs.
        :param str pipeline_type: project or sample-level pipeline
        :return dict: {table_name: model}
        """
        _LOGGER.debug(f"Creating models for '{self.pipeline_name}' table in '{PKG_NAME}' database")
        model = self.parsed_schema.build_model(pipeline_type=pipeline_type)
        table_name = self.parsed_schema._table_name(pipeline_type)
        # TODO reconsider line below. Why do we need to return a dict?
        if model:
            return {table_name: model}
        else:
            raise SchemaError(
                f"Neither project nor samples model could be built from schema source: {self.status_schema_source}"
            )

    def _create_history_orms(self, pipeline_type):
        """Creates the additional ORMs for auditing result modifications
        :param str pipeline_type: project or sample-level pipeline
        :return dict: {table_name: model}
        """
        model, table_name = self.parsed_schema.build_history_model(pipeline_type=pipeline_type)
        if model:
            return {table_name: model}
        else:
            raise SchemaError(
                f"Neither project nor samples model could be built from schema source: {self.status_schema_source}"
            )

    @property
    def _engine(self):
        """Access the database engine backing this manager."""
        try:
            return self.db_engine_key
        except AttributeError:
            # Do it this way rather than .setdefault to avoid evaluating
            # the expression for the default argument (i.e., building
            # the engine) if it's not necessary.
            self.db_engine_key = create_engine(self.db_url, echo=self.show_db_logs)
            return self.db_engine_key

    @property
    @contextmanager
    def session(self):
        """
        Provide a transactional scope around a series of query
        operations.
        """
        session = Session(self._engine)
        _LOGGER.debug("Created session")
        try:
            yield session
        except:
            _LOGGER.info("session.rollback")
            session.rollback()
            raise
        finally:
            # _LOGGER.info("session.close")
            session.close()
        _LOGGER.debug("Ending session")
