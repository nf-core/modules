from typing import Optional


class BasePephubclientException(Exception):
    def __init__(self, message: str):
        super().__init__(message)


class IncorrectQueryStringError(BasePephubclientException):
    def __init__(self, query_string: Optional[str] = None):
        self.query_string = query_string
        super().__init__(
            f"PEP data with passed namespace and project ({self.query_string}) name not found."
        )


class ResponseError(BasePephubclientException):
    default_message = "The response looks incorrect and must be verified manually."

    def __init__(self, message: Optional[str] = None):
        self.message = message
        super().__init__(self.message or self.default_message)


class PEPExistsError(BasePephubclientException):
    default_message = (
        "PEP already exists. Change location, delete previous PEP, or set force argument "
        "to overwrite previous PEP"
    )

    def __init__(self, message: Optional[str] = None):
        self.message = message
        super().__init__(self.message or self.default_message)
