from enum import Enum
from typing import Optional
import os

from pydantic import BaseModel, field_validator

PEPHUB_BASE_URL = os.getenv(
    "PEPHUB_BASE_URL", default="https://pephub-api.databio.org/"
)
# PEPHUB_BASE_URL = "http://0.0.0.0:8000/"
PEPHUB_PEP_API_BASE_URL = f"{PEPHUB_BASE_URL}api/v1/projects/"
PEPHUB_PEP_SEARCH_URL = f"{PEPHUB_BASE_URL}api/v1/namespaces/{{namespace}}/projects"
PEPHUB_PUSH_URL = f"{PEPHUB_BASE_URL}api/v1/namespaces/{{namespace}}/projects/json"

PEPHUB_SAMPLE_URL = f"{PEPHUB_BASE_URL}api/v1/projects/{{namespace}}/{{project}}/samples/{{sample_name}}"
PEPHUB_VIEW_URL = (
    f"{PEPHUB_BASE_URL}api/v1/projects/{{namespace}}/{{project}}/views/{{view_name}}"
)
PEPHUB_VIEW_SAMPLE_URL = f"{PEPHUB_BASE_URL}api/v1/projects/{{namespace}}/{{project}}/views/{{view_name}}/{{sample_name}}"


class RegistryPath(BaseModel):
    protocol: Optional[str] = None
    namespace: str
    item: str
    subitem: Optional[str] = None
    tag: Optional[str] = "default"

    @field_validator("tag")
    def tag_should_not_be_none(cls, v):
        return v or "default"


class ResponseStatusCodes(int, Enum):
    OK = 200
    ACCEPTED = 202
    UNAUTHORIZED = 401
    FORBIDDEN = 403
    NOT_EXIST = 404
    CONFLICT = 409
    INTERNAL_ERROR = 500


USER_DATA_FILE_NAME = "jwt.txt"
HOME_PATH = os.getenv("HOME")
if not HOME_PATH:
    HOME_PATH = os.path.expanduser("~")
PATH_TO_FILE_WITH_JWT = os.path.join(HOME_PATH, ".pephubclient/") + USER_DATA_FILE_NAME
