import datetime
from typing import Optional, List, Union

from pydantic import BaseModel, Field, field_validator, ConfigDict
from peppy.const import CONFIG_KEY, SUBSAMPLE_RAW_LIST_KEY, SAMPLE_RAW_DICT_KEY


class ProjectDict(BaseModel):
    """
    Project dict (raw) model
    """

    config: dict = Field(alias=CONFIG_KEY)
    subsample_list: Optional[list] = Field(alias=SUBSAMPLE_RAW_LIST_KEY)
    sample_list: list = Field(alias=SAMPLE_RAW_DICT_KEY)

    model_config = ConfigDict(populate_by_name=True, extra="allow")


class ProjectUploadData(BaseModel):
    """
    Model used in post request to upload project
    """

    pep_dict: ProjectDict
    tag: Optional[str] = "default"
    is_private: Optional[bool] = False
    overwrite: Optional[bool] = False

    @field_validator("tag")
    def tag_should_not_be_none(cls, v):
        return v or "default"


class ProjectAnnotationModel(BaseModel):
    namespace: str
    name: str
    tag: str
    is_private: bool
    number_of_samples: int
    description: str
    last_update_date: datetime.datetime
    submission_date: datetime.datetime
    digest: str
    pep_schema: Union[str, int, None] = None
    pop: bool = False
    stars_number: Optional[int] = 0
    forked_from: Optional[Union[str, None]] = None


class SearchReturnModel(BaseModel):
    count: int
    limit: int
    offset: int
    results: List[ProjectAnnotationModel]
