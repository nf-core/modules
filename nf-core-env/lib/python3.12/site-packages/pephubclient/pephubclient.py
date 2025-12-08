from typing import NoReturn, Optional, Literal
from typing_extensions import deprecated

import peppy
from peppy.const import NAME_KEY
import urllib3
from pydantic import ValidationError
from ubiquerg import parse_registry_path

from pephubclient.constants import (
    PEPHUB_PEP_API_BASE_URL,
    PEPHUB_PUSH_URL,
    RegistryPath,
    ResponseStatusCodes,
    PEPHUB_PEP_SEARCH_URL,
    PATH_TO_FILE_WITH_JWT,
)
from pephubclient.exceptions import (
    IncorrectQueryStringError,
    ResponseError,
)
from pephubclient.files_manager import FilesManager
from pephubclient.helpers import MessageHandler, RequestManager, save_pep
from pephubclient.models import (
    ProjectDict,
    ProjectUploadData,
    SearchReturnModel,
    ProjectAnnotationModel,
)
from pephubclient.pephub_oauth.pephub_oauth import PEPHubAuth
from pephubclient.modules.view import PEPHubView
from pephubclient.modules.sample import PEPHubSample

urllib3.disable_warnings()


class PEPHubClient(RequestManager):
    def __init__(self):
        self.__jwt_data = FilesManager.load_jwt_data_from_file(PATH_TO_FILE_WITH_JWT)

        self.__view = PEPHubView(self.__jwt_data)
        self.__sample = PEPHubSample(self.__jwt_data)

    @property
    def view(self) -> PEPHubView:
        return self.__view

    @property
    def sample(self) -> PEPHubSample:
        return self.__sample

    def login(self) -> NoReturn:
        """
        Log in to PEPhub
        """
        user_token = PEPHubAuth().login_to_pephub()

        FilesManager.save_jwt_data_to_file(PATH_TO_FILE_WITH_JWT, user_token)
        self.__jwt_data = FilesManager.load_jwt_data_from_file(PATH_TO_FILE_WITH_JWT)

    def logout(self) -> NoReturn:
        """
        Log out from PEPhub
        """
        FilesManager.delete_file_if_exists(PATH_TO_FILE_WITH_JWT)
        self.__jwt_data = None

    def pull(
        self,
        project_registry_path: str,
        force: Optional[bool] = False,
        zip: Optional[bool] = False,
        output: Optional[str] = None,
    ) -> None:
        """
        Download project locally

        :param str project_registry_path: Project registry path in PEPhub (e.g. databio/base:default)
        :param bool force: if project exists, overwrite it.
        :param bool zip: if True, save project as zip file
        :param str output: path where project will be saved
        :return: None
        """
        project_dict = self.load_raw_pep(
            registry_path=project_registry_path,
        )

        save_pep(
            project=project_dict,
            reg_path=project_registry_path,
            force=force,
            project_path=output,
            zip=zip,
        )

    def load_project(
        self,
        project_registry_path: str,
        query_param: Optional[dict] = None,
    ) -> peppy.Project:
        """
        Load peppy project from PEPhub in peppy.Project object

        :param project_registry_path: registry path of the project
        :param query_param: query parameters used in get request
        :return Project: peppy project.
        """
        raw_pep = self.load_raw_pep(project_registry_path, query_param)
        peppy_project = peppy.Project().from_dict(raw_pep)
        return peppy_project

    def push(
        self,
        cfg: str,
        namespace: str,
        name: Optional[str] = None,
        tag: Optional[str] = None,
        is_private: Optional[bool] = False,
        force: Optional[bool] = False,
    ) -> None:
        """
        Push (upload/update) project to Pephub using config/csv path

        :param str cfg: Project config file (YAML) or sample table (CSV/TSV)
            with one row per sample to constitute project
        :param str namespace: namespace
        :param str name: project name
        :param str tag: project tag
        :param bool is_private: Specifies whether project should be private [Default= False]
        :param bool force: Force push to the database. Use it to update, or upload project. [Default= False]
        :return: None
        """
        peppy_project = peppy.Project(cfg=cfg)
        self.upload(
            project=peppy_project,
            namespace=namespace,
            name=name,
            tag=tag,
            is_private=is_private,
            force=force,
        )

    def upload(
        self,
        project: peppy.Project,
        namespace: str,
        name: str = None,
        tag: str = None,
        is_private: bool = False,
        force: bool = True,
    ) -> None:
        """
        Upload peppy project to the PEPhub.

        :param peppy.Project project: Project object that has to be uploaded to the DB
        :param namespace: namespace
        :param name: project name
        :param tag: project tag
        :param force: Force push to the database. Use it to update, or upload project.
        :param is_private: Make project private
        :param force: overwrite project if it exists
        :return: None
        """
        if name:
            project[NAME_KEY] = name

        upload_data = ProjectUploadData(
            pep_dict=project.to_dict(
                extended=True,
                orient="records",
            ),
            tag=tag,
            is_private=is_private,
            overwrite=force,
        )
        pephub_response = self.send_request(
            method="POST",
            url=self._build_push_request_url(namespace=namespace),
            headers=self.parse_header(self.__jwt_data),
            json=upload_data.model_dump(),
            cookies=None,
        )
        if pephub_response.status_code == ResponseStatusCodes.ACCEPTED:
            MessageHandler.print_success(
                f"Project '{namespace}/{name}:{upload_data.tag}' was successfully uploaded"
            )
        elif pephub_response.status_code == ResponseStatusCodes.CONFLICT:
            raise ResponseError(
                "Project already exists. Set force to overwrite project."
            )
        elif pephub_response.status_code == ResponseStatusCodes.UNAUTHORIZED:
            raise ResponseError("Unauthorized! Failure in uploading project.")
        elif pephub_response.status_code == ResponseStatusCodes.FORBIDDEN:
            raise ResponseError(
                "User does not have permission to write to this namespace!"
            )
        else:
            raise ResponseError(
                f"Unexpected Response Error. {pephub_response.status_code}"
            )
        return None

    def find_project(
        self,
        namespace: str,
        query_string: str = "",
        limit: int = 100,
        offset: int = 0,
        filter_by: Literal["submission_date", "last_update_date"] = None,
        start_date: str = None,
        end_date: str = None,
    ) -> SearchReturnModel:
        """
        Find project in specific namespace and return list of PEP annotation

        :param namespace: Namespace where to search for projects
        :param query_string: Search query
        :param limit: Return limit
        :param offset: Return offset
        :param filter_by: Use filter date. Option: [submission_date, last_update_date]
        :param start_date: filter beginning date
        :param end_date: filter end date (if none today's date is used)
        :return:
        """

        query_param = {
            "q": query_string,
            "limit": limit,
            "offset": offset,
        }
        if filter_by in ["submission_date", "last_update_date"]:
            query_param["filter_by"] = filter_by
            query_param["filter_start_date"] = start_date
            if end_date:
                query_param["filter_end_date"] = end_date

        url = self._build_project_search_url(
            namespace=namespace,
            query_param=query_param,
        )

        pephub_response = self.send_request(
            method="GET",
            url=url,
            headers=self.parse_header(self.__jwt_data),
            json=None,
            cookies=None,
        )
        if pephub_response.status_code == ResponseStatusCodes.OK:
            decoded_response = self.decode_response(pephub_response, output_json=True)
            project_list = []
            for project_found in decoded_response["results"]:
                project_list.append(ProjectAnnotationModel(**project_found))
            return SearchReturnModel(**decoded_response)

    @deprecated("This method is deprecated. Use load_raw_pep instead.")
    def _load_raw_pep(
        self,
        registry_path: str,
        jwt_data: Optional[str] = None,
        query_param: Optional[dict] = None,
    ) -> dict:
        """
        !!! This method is deprecated. Use load_raw_pep instead. !!!

        Request PEPhub and return the requested project as peppy.Project object.

        :param registry_path: Project namespace, eg. "geo/GSE124224:tag"
        :param query_param: Optional variables to be passed to PEPhub
        :return: Raw project in dict.
        """
        return self.load_raw_pep(registry_path, query_param)

    def load_raw_pep(
        self,
        registry_path: str,
        query_param: Optional[dict] = None,
    ) -> dict:
        """
        Request PEPhub and return the requested project as peppy.Project object.

        :param registry_path: Project namespace, eg. "geo/GSE124224:tag"
        :param query_param: Optional variables to be passed to PEPhub
        :return: Raw project in dict.
        """
        query_param = query_param or {}
        query_param["raw"] = "true"

        self._set_registry_data(registry_path)
        pephub_response = self.send_request(
            method="GET",
            url=self._build_pull_request_url(query_param=query_param),
            headers=self.parse_header(self.__jwt_data),
            cookies=None,
        )
        if pephub_response.status_code == ResponseStatusCodes.OK:
            decoded_response = self.decode_response(pephub_response, output_json=True)
            correct_proj_dict = ProjectDict(**decoded_response)

            # This step is necessary because of this issue: https://github.com/pepkit/pephub/issues/124
            return correct_proj_dict.model_dump(by_alias=True)

        if pephub_response.status_code == ResponseStatusCodes.NOT_EXIST:
            raise ResponseError("File does not exist, or you are unauthorized.")
        if pephub_response.status_code == ResponseStatusCodes.INTERNAL_ERROR:
            raise ResponseError(
                f"Internal server error. Unexpected return value. Error: {pephub_response.status_code}"
            )

    def _set_registry_data(self, query_string: str) -> None:
        """
        Parse provided query string to extract project name, sample name, etc.

        :param query_string: Passed by user. Contain information needed to locate the project.
        :return: Parsed query string.
        """
        try:
            self.registry_path = RegistryPath(**parse_registry_path(query_string))
        except (ValidationError, TypeError):
            raise IncorrectQueryStringError(query_string=query_string)

    def _build_pull_request_url(self, query_param: dict = None) -> str:
        """
        Build request for getting projects form pephub

        :param query_param: dict of parameters used in query string
        :return: url string
        """
        query_param = query_param or {}
        query_param["tag"] = self.registry_path.tag

        endpoint = self.registry_path.namespace + "/" + self.registry_path.item

        variables_string = self.parse_query_param(query_param)
        endpoint += variables_string

        return PEPHUB_PEP_API_BASE_URL + endpoint

    @staticmethod
    def _build_project_search_url(namespace: str, query_param: dict = None) -> str:
        """
        Build request for searching projects form pephub

        :param query_param: dict of parameters used in query string
        :return: url string
        """

        variables_string = RequestManager.parse_query_param(query_param)
        endpoint = variables_string

        return PEPHUB_PEP_SEARCH_URL.format(namespace=namespace) + endpoint

    @staticmethod
    def _build_push_request_url(namespace: str) -> str:
        """
        Build project uplaod request used in pephub

        :param namespace: namespace where project will be uploaded
        :return: url string
        """
        return PEPHUB_PUSH_URL.format(namespace=namespace)
