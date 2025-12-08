import logging

from pephubclient.helpers import RequestManager
from pephubclient.constants import PEPHUB_SAMPLE_URL, ResponseStatusCodes
from pephubclient.exceptions import ResponseError

_LOGGER = logging.getLogger("pephubclient")


class PEPHubSample(RequestManager):
    """
    Class for managing samples in PEPhub and provides methods for
        getting, creating, updating and removing samples.
    This class is not related to peppy.Sample class.
    """

    def __init__(self, jwt_data: str = None):
        """
        :param jwt_data: jwt token for authorization
        """

        self.__jwt_data = jwt_data

    def get(
        self,
        namespace: str,
        name: str,
        tag: str,
        sample_name: str = None,
    ) -> dict:
        """
        Get sample from project in PEPhub.

        :param namespace: namespace of project
        :param name: name of project
        :param tag: tag of project
        :param sample_name: sample name
        :return: Sample object
        """
        url = self._build_sample_request_url(
            namespace=namespace, name=name, sample_name=sample_name
        )

        url = url + self.parse_query_param(pep_variables={"tag": tag})

        response = self.send_request(
            method="GET", url=url, headers=self.parse_header(self.__jwt_data)
        )
        if response.status_code == ResponseStatusCodes.OK:
            return self.decode_response(response, output_json=True)
        if response.status_code == ResponseStatusCodes.NOT_EXIST:
            raise ResponseError(
                f"Sample does not exist. Project: '{namespace}/{name}:{tag}'. Sample_name: '{sample_name}'"
            )
        elif response.status_code == ResponseStatusCodes.INTERNAL_ERROR:
            raise ResponseError("Internal server error. Unexpected return value.")
        else:
            raise ResponseError(
                f"Unexpected return value. Error: {response.status_code}"
            )

    def create(
        self,
        namespace: str,
        name: str,
        tag: str,
        sample_name: str,
        sample_dict: dict,
        overwrite: bool = False,
    ) -> None:
        """
        Create sample in project in PEPhub.

        :param namespace: namespace of project
        :param name: name of project
        :param tag: tag of project
        :param sample_dict: sample dict
        :param sample_name: sample name
        :param overwrite: overwrite sample if it exists
        :return: None
        """
        url = self._build_sample_request_url(
            namespace=namespace,
            name=name,
            sample_name=sample_name,
        )

        url = url + self.parse_query_param(
            pep_variables={"tag": tag, "overwrite": overwrite}
        )

        # add sample name to sample_dict if it is not there
        if sample_name not in sample_dict.values():
            sample_dict["sample_name"] = sample_name

        response = self.send_request(
            method="POST",
            url=url,
            headers=self.parse_header(self.__jwt_data),
            json=sample_dict,
        )
        if response.status_code == ResponseStatusCodes.ACCEPTED:
            _LOGGER.info(
                f"Sample '{sample_name}' added to project '{namespace}/{name}:{tag}' successfully."
            )
            return None
        elif response.status_code == ResponseStatusCodes.NOT_EXIST:
            raise ResponseError(f"Project '{namespace}/{name}:{tag}' does not exist.")
        elif response.status_code == ResponseStatusCodes.CONFLICT:
            raise ResponseError(
                f"Sample '{sample_name}' already exists. Set overwrite to True to overwrite sample."
            )
        else:
            raise ResponseError(
                f"Unexpected return value. Error: {response.status_code}"
            )

    def update(
        self,
        namespace: str,
        name: str,
        tag: str,
        sample_name: str,
        sample_dict: dict,
    ):
        """
        Update sample in project in PEPhub.

        :param namespace: namespace of project
        :param name: name of project
        :param tag: tag of project
        :param sample_name: sample name
        :param sample_dict: sample dict, that contain elements to update, or
        :return: None
        """

        url = self._build_sample_request_url(
            namespace=namespace, name=name, sample_name=sample_name
        )

        url = url + self.parse_query_param(pep_variables={"tag": tag})

        response = self.send_request(
            method="PATCH",
            url=url,
            headers=self.parse_header(self.__jwt_data),
            json=sample_dict,
        )
        if response.status_code == ResponseStatusCodes.ACCEPTED:
            _LOGGER.info(
                f"Sample '{sample_name}' updated in project '{namespace}/{name}:{tag}' successfully."
            )
            return None
        elif response.status_code == ResponseStatusCodes.NOT_EXIST:
            raise ResponseError(
                f"Sample '{sample_name}' or project {namespace}/{name}:{tag} does not exist. Error: {response.status_code}"
            )
        else:
            raise ResponseError(
                f"Unexpected return value. Error: {response.status_code}"
            )

    def remove(self, namespace: str, name: str, tag: str, sample_name: str):
        """
        Remove sample from project in PEPhub.

        :param namespace: namespace of project
        :param name: name of project
        :param tag: tag of project
        :param sample_name: sample name
        :return: None
        """
        url = self._build_sample_request_url(
            namespace=namespace, name=name, sample_name=sample_name
        )

        url = url + self.parse_query_param(pep_variables={"tag": tag})

        response = self.send_request(
            method="DELETE",
            url=url,
            headers=self.parse_header(self.__jwt_data),
        )
        if response.status_code == ResponseStatusCodes.ACCEPTED:
            _LOGGER.info(
                f"Sample '{sample_name}' removed from project '{namespace}/{name}:{tag}' successfully."
            )
            return None
        elif response.status_code == ResponseStatusCodes.NOT_EXIST:
            raise ResponseError(
                f"Sample '{sample_name}' or project {namespace}/{name}:{tag} does not exist. Error: {response.status_code}"
            )
        else:
            raise ResponseError(
                f"Unexpected return value. Error: {response.status_code}"
            )

    @staticmethod
    def _build_sample_request_url(namespace: str, name: str, sample_name: str) -> str:
        """
        Build url for sample request.

        :param namespace: namespace where project will be uploaded
        :return: url string
        """
        return PEPHUB_SAMPLE_URL.format(
            namespace=namespace, project=name, sample_name=sample_name
        )
