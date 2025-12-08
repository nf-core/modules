import json
from typing import Any, Callable, Optional, Union
import peppy
import yaml
import os
import pandas as pd
from peppy.const import (
    NAME_KEY,
    DESC_KEY,
    CONFIG_KEY,
    SUBSAMPLE_RAW_LIST_KEY,
    SAMPLE_RAW_DICT_KEY,
    CFG_SAMPLE_TABLE_KEY,
    CFG_SUBSAMPLE_TABLE_KEY,
)

import requests
from requests.exceptions import ConnectionError
from urllib.parse import urlencode

from ubiquerg import parse_registry_path
from pydantic import ValidationError

from pephubclient.exceptions import PEPExistsError, ResponseError
from pephubclient.constants import RegistryPath
from pephubclient.files_manager import FilesManager
from pephubclient.models import ProjectDict


class RequestManager:
    @staticmethod
    def send_request(
        method: str,
        url: str,
        headers: Optional[dict] = None,
        cookies: Optional[dict] = None,
        params: Optional[dict] = None,
        json: Optional[Union[dict, list]] = None,
    ) -> requests.Response:
        request_return = requests.request(
            method=method,
            url=url,
            verify=False,
            cookies=cookies,
            headers=headers,
            params=params,
            json=json,
            timeout=10,
        )
        if request_return.status_code == 401:
            if (
                RequestManager.decode_response(request_return, output_json=True).get(
                    "detail"
                )
                == "JWT has expired"
            ):
                raise ResponseError("JWT has expired. Please log in again.")
        return request_return

    @staticmethod
    def decode_response(
        response: requests.Response, encoding: str = "utf-8", output_json: bool = False
    ) -> Union[str, dict]:
        """
        Decode the response from GitHub and pack the returned data into appropriate model.

        :param response: Response from GitHub.
        :param encoding: Response encoding [Default: utf-8]
        :param output_json: If True, return response in json format
        :return: Response data as an instance of correct model.
        """

        try:
            if output_json:
                return response.json()
            else:
                return response.content.decode(encoding)
        except json.JSONDecodeError as err:
            raise ResponseError(f"Error in response encoding format: {err}")

    @staticmethod
    def parse_query_param(pep_variables: dict) -> str:
        """
        Grab all the variables passed by user (if any) and parse them to match the format specified
        by PEPhub API for query parameters.

        :param pep_variables: dict of query parameters
        :return: PEPHubClient variables transformed into string in correct format.
        """
        return "?" + urlencode(pep_variables)

    @staticmethod
    def parse_header(jwt_data: Optional[str] = None) -> dict:
        """
        Create Authorization header

        :param jwt_data: jwt string
        :return: Authorization dict
        """
        if jwt_data:
            return {"Authorization": jwt_data}
        else:
            return {}


class MessageHandler:
    """
    Class holding print function in different colors
    """

    RED = 9
    YELLOW = 11
    GREEN = 40

    @staticmethod
    def print_error(text: str) -> None:
        print(f"\033[38;5;9m{text}\033[0m")

    @staticmethod
    def print_success(text: str) -> None:
        print(f"\033[38;5;40m{text}\033[0m")

    @staticmethod
    def print_warning(text: str) -> None:
        print(f"\033[38;5;11m{text}\033[0m")


def call_client_func(func: Callable[..., Any], **kwargs) -> Any:
    """
    Catch exceptions in functions called through cli.

    :param func: The function to call.
    :param kwargs: The keyword arguments to pass to the function.
    :return: The result of the function call.
    """

    try:
        func(**kwargs)
    except ConnectionError as err:
        MessageHandler.print_error(f"Failed to connect to server. Try later. {err}")
    except ResponseError as err:
        MessageHandler.print_error(f"{err}")
    except PEPExistsError as err:
        MessageHandler.print_warning(f"PEP already exists. {err}")
    except OSError as err:
        MessageHandler.print_error(f"{err}")


def is_registry_path(input_string: str) -> bool:
    """
    Check if input is a registry path to pephub
    :param str input_string: path to the PEP (or registry path)
    :return bool: True if input is a registry path
    """
    if input_string.endswith(".yaml"):
        return False
    try:
        RegistryPath(**parse_registry_path(input_string))
    except (ValidationError, TypeError):
        return False
    return True


def unwrap_registry_path(input_string: str) -> RegistryPath:
    """
    Unwrap registry path from string
    :param str input_string: path to the PEP (or registry path)
    :return RegistryPath: RegistryPath object
    """
    return RegistryPath(**parse_registry_path(input_string))


def _build_filename(registry_path: RegistryPath) -> str:
    """
    Takes query string and creates output filename to save the project to.

    :param registry_path: Query string that was used to find the project.
    :return: Filename uniquely identifying the project.
    """
    filename = "_".join(filter(bool, [registry_path.namespace, registry_path.item]))
    if registry_path.tag:
        filename += f"_{registry_path.tag}"
    return filename


def _save_zip_pep(project: dict, zip_filepath: str, force: bool = False) -> None:
    """
    Zip and save a project

    :param project: peppy project to zip
    :param zip_filepath: path to save zip file
    :param force: overwrite project if exists
    """

    content_to_zip = {}
    config = project[CONFIG_KEY]
    project_name = config[NAME_KEY]

    if project[SAMPLE_RAW_DICT_KEY] is not None:
        config[CFG_SAMPLE_TABLE_KEY] = ["sample_table.csv"]
        content_to_zip["sample_table.csv"] = pd.DataFrame(
            project[SAMPLE_RAW_DICT_KEY]
        ).to_csv(index=False)

    if project[SUBSAMPLE_RAW_LIST_KEY] is not None:
        if not isinstance(project[SUBSAMPLE_RAW_LIST_KEY], list):
            config[CFG_SUBSAMPLE_TABLE_KEY] = ["subsample_table1.csv"]
            content_to_zip["subsample_table1.csv"] = pd.DataFrame(
                project[SUBSAMPLE_RAW_LIST_KEY]
            ).to_csv(index=False)
        else:
            config[CFG_SUBSAMPLE_TABLE_KEY] = []
            for number, file in enumerate(project[SUBSAMPLE_RAW_LIST_KEY]):
                file_name = f"subsample_table{number + 1}.csv"
                config[CFG_SUBSAMPLE_TABLE_KEY].append(file_name)
                content_to_zip[file_name] = pd.DataFrame(file).to_csv(index=False)

    content_to_zip[f"{project_name}_config.yaml"] = yaml.dump(config, indent=4)
    FilesManager.save_zip_file(content_to_zip, file_path=zip_filepath, force=force)

    MessageHandler.print_success(f"Project was saved successfully -> {zip_filepath}")
    return None


def _save_unzipped_pep(
    project_dict: dict, folder_path: str, force: bool = False
) -> None:
    """
    Save unzipped project to specified folder

    :param project_dict: raw pep project
    :param folder_path: path to save project
    :param force: overwrite project if exists
    :return: None
    """

    def full_path(fn: str) -> str:
        return os.path.join(folder_path, fn)

    project_name = project_dict[CONFIG_KEY][NAME_KEY]
    sample_table_filename = "sample_table.csv"
    yaml_full_path = full_path(f"{project_name}_config.yaml")
    sample_full_path = full_path(sample_table_filename)
    if not force:
        extant = [p for p in [yaml_full_path, sample_full_path] if os.path.isfile(p)]
        if extant:
            raise PEPExistsError(f"{len(extant)} file(s) exist(s): {', '.join(extant)}")

    config_dict = project_dict.get(CONFIG_KEY)
    config_dict[NAME_KEY] = project_name
    config_dict[DESC_KEY] = project_dict[CONFIG_KEY][DESC_KEY]
    config_dict["sample_table"] = sample_table_filename

    sample_pandas = pd.DataFrame(project_dict.get(SAMPLE_RAW_DICT_KEY, {}))

    subsample_list = [
        pd.DataFrame(sub_a) for sub_a in project_dict.get(SUBSAMPLE_RAW_LIST_KEY) or []
    ]

    filenames = []
    for idx, subsample in enumerate(subsample_list):
        fn = f"subsample_table{idx + 1}.csv"
        filenames.append(fn)
        FilesManager.save_pandas(subsample, full_path(fn), not_force=False)
    config_dict["subsample_table"] = filenames

    FilesManager.save_yaml(config_dict, yaml_full_path, not_force=False)
    FilesManager.save_pandas(sample_pandas, sample_full_path, not_force=False)

    if config_dict.get("subsample_table"):
        for number, subsample in enumerate(subsample_list):
            FilesManager.save_pandas(
                subsample,
                os.path.join(folder_path, config_dict["subsample_table"][number]),
                not_force=False,
            )

    MessageHandler.print_success(f"Project was saved successfully -> {folder_path}")
    return None


def save_pep(
    project: Union[dict, peppy.Project],
    reg_path: str = None,
    force: bool = False,
    project_path: Optional[str] = None,
    zip: bool = False,
) -> None:
    """
    Save project locally.

    :param dict project: PEP dictionary (raw project)
    :param str reg_path: Project registry path in PEPhub (e.g. databio/base:default). If not provided,
        folder will be created with just project name.
    :param bool force: overwrite project if exists
    :param str project_path: Path where project will be saved. By default, it will be saved in current directory.
    :param bool zip: If True, save project as zip file
    :return: None
    """
    if isinstance(project, peppy.Project):
        project = project.to_dict(extended=True, orient="records")

    project = ProjectDict(**project).model_dump(by_alias=True)

    if not project_path:
        project_path = os.getcwd()

    if reg_path:
        file_name = _build_filename(RegistryPath(**parse_registry_path(reg_path)))
    else:
        file_name = project[CONFIG_KEY][NAME_KEY]

    if zip:
        _save_zip_pep(
            project,
            zip_filepath=f"{os.path.join(project_path, file_name)}.zip",
            force=force,
        )
        return None

    folder_path = FilesManager.create_project_folder(
        parent_path=project_path, folder_name=file_name
    )
    _save_unzipped_pep(project, folder_path, force=force)
