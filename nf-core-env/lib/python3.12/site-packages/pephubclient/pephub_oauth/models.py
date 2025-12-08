from pydantic import BaseModel


class InitializeDeviceCodeResponse(BaseModel):
    device_code: str
    auth_url: str


class PEPHubDeviceTokenResponse(BaseModel):
    jwt_token: str
