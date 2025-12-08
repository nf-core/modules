import logging
import os
from typing import List, Optional, Tuple, Union

import fastapi
import uvicorn
from pydantic import BaseModel

from pipestat import SamplePipestatManager
from pipestat.exceptions import RecordNotFoundError
from pipestat.reports import fetch_pipeline_results

_LOGGER = logging.getLogger(__name__)

app = fastapi.FastAPI(
    title="Pipestat Reader",
    description="Allows reading pipestat files from a PostgreSQL Database",
    version="0.01",
)


class FilterQuery(BaseModel):
    column_names: Union[List[str], None] = None
    filter_conditions: Union[List[Tuple[str, str, str]], None] = None


@app.get("/")
async def home():
    """
    Display the home page
    """
    return {"Home": "Welcome"}


@app.get("/data/{record_identifier}")
async def retrieve_results_one_record(record_identifier: str):
    """
    Get all the results for one record
    """
    try:
        result = psm.retrieve(record_identifier=record_identifier)
    except RecordNotFoundError:
        return {"result": "Record not found"}
    return {"result": result}


@app.get("/data/")
async def retrieve_all_records():
    """
    Get all reported records
    """
    try:
        result = psm.get_records()
    except RecordNotFoundError:
        return {"result": "Record not found"}
    return {"result": result}


@app.get("/data/{record_identifier}/{result_identifier}")
async def retrieve_results(record_identifier: str, result_identifier: str):
    """
    Get specific result given a record identifier and a result identifier
    """
    try:
        result = psm.retrieve(
            record_identifier=record_identifier, result_identifier=result_identifier
        )
    except RecordNotFoundError:
        # TODO this should be more specific than record not found because maybe it's just the result that does not exist
        return {"result": "Record not found"}
    return {"result": result}


@app.get("/output_schema/")
async def retrieve_output_schema(pipeline_type: Optional[str] = None):
    """
    Get the output_schema used by the PipestatManager.
    """

    if pipeline_type == "sample":
        return {"output schema": psm.schema._sample_level_data}
    if pipeline_type == "project":
        return {"output schema": psm.schema._project_level_data}
    if pipeline_type is None:
        return {"output schema": psm.schema}
    else:
        return {"output schema": "output schema not found"}


@app.get("/all_table_contents/")
async def retrieve_table_contents():
    """
    Get all table contents
    """
    # Add skip and limit here as well.
    results = psm.backend.select()
    return {"table_contents": results}


@app.get("/{file_type}/")
async def retrieve_filetype(file_type: str):
    """
    Get all records by filetype
    """
    records_by_filetype = []
    all_records = psm.get_records()["records"]
    for sample in all_records:
        file_result = fetch_pipeline_results(
            project=psm,
            sample_name=sample,
            inclusion_fun=lambda x: x == file_type,
        )
        records_by_filetype.append(file_result)
    return {"records_by_filetype": records_by_filetype}


@app.post("/filtered_table_contents/")
async def retrieve_filtered_table_contents(query_filter: Optional[FilterQuery] = None):
    """
        Get column contents for specific column names and/or filter conditions

    {
      "column_names": [
        "md5sum", "status"
      ],
      "filter_conditions": [["record_identifier", "eq", "random_sample_id2"]]
    }

    """
    try:
        results = psm.backend.select(
            columns=query_filter.column_names,
            filter_conditions=query_filter.filter_conditions,
        )
    except AttributeError:
        return {"response": f"Attribute error for query: {query_filter.column_names}"}
    return {"response": results}


def create_global_pipestatmanager(pipestatcfg):
    """
    build a global pipestatmanager to be used by the endpoints
    """
    global psm
    psm = SamplePipestatManager(config_file=pipestatcfg)
    print("GLOBAL PSM CREATED")


if __name__ != "__main__":
    if os.environ.get("PIPESTAT_CONFIG") is not None:
        pipestatcfg = os.environ.get("PIPESTAT_CONFIG")
        create_global_pipestatmanager(pipestatcfg)
    else:
        _LOGGER.error("Configure by setting PIPESTAT_CONFIG env var")


def main(
    configfile: Optional[str] = None,
    host: Optional[str] = None,
    port: Optional[int] = None,
):
    """
    passes relevant info to create a global pipestat manager and then utilizes uvicorn to run at the host address and
    port. These parameters are passed during `pipestat serve` cli usage.
    param: str configfile: a path to a pipestat configfile
    param: str host: host address for uvicorn server
    param: str port: port number for uvicorn server

    """
    pipestatcfg = configfile or os.environ.get("PIPESTAT_CONFIG")
    if pipestatcfg is None:
        _LOGGER.error(
            "A pipestat configuration file must be supplied."
            "Configure by setting PIPESTAT_CONFIG env var or pass as input argument using --config"
        )
    create_global_pipestatmanager(pipestatcfg)
    # Note input argument app vs "reader:app" causes different behavior when using uvicorn.run
    uvicorn.run(
        app,
        host=host or "0.0.0.0",
        port=port or 8001,
    )


if __name__ == "__main__":
    main()
