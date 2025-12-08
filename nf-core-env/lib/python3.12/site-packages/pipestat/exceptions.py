"""Package exception types"""

from typing import Iterable, Optional

from .const import CLASSES_BY_TYPE, ENV_VARS

__all__ = [
    "InvalidTypeError",
    "IncompatibleClassError",
    "NoBackendSpecifiedError",
    "PipestatError",
    "PipestatDatabaseError",
    "MissingConfigDataError",
    "SchemaError",
    "SchemaNotFoundError",
    "PipestatDataError",
    "UnrecognizedStatusError",
    "RecordNotFoundError",
    "PipelineTypeNotSuppliedError",
    "InvalidTimeFormatError",
    "PipestatDependencyError",
    "ColumnNotFoundError",
    "SchemaValidationErrorDuringReport",
    "PipestatPEPHubError",
    "PipestatSummarizeError",
]


class RecordNotFoundError(LookupError):
    def __init__(self, msg):
        super(RecordNotFoundError, self).__init__(msg)


class ColumnNotFoundError(LookupError):
    """A specified attribute (column) is not in the table or schema"""

    def __init__(self, msg):
        super(ColumnNotFoundError, self).__init__(msg)


class PipelineTypeNotSuppliedError(LookupError):
    def __init__(self, msg):
        super(PipelineTypeNotSuppliedError, self).__init__(msg)


class PipestatError(Exception):
    """Base exception type for this package"""


class NoBackendSpecifiedError(PipestatError):
    """Subtype for designating lack of backend specification"""


class SchemaError(PipestatError):
    """Schema error"""

    def __init__(self, msg):
        super(SchemaError, self).__init__(msg)


class SchemaNotFoundError(SchemaError):
    """Schema not found error"""

    def __init__(self, msg, cli=False):
        txt = f"Results schema not found. The schema is required to {msg}. "
        txt += (
            "It needs to be supplied as an CLI argument"
            if cli
            else "It needs to be supplied to the object constructor"
        )
        txt += f" or via '{ENV_VARS['schema']}' environment variable."
        super(SchemaNotFoundError, self).__init__(txt)


class SchemaValidationErrorDuringReport(SchemaError):
    """Adds clarity to JSON schema validation errors by providing additional information to error message."""

    def __init__(self, msg, record_identifier, result_identifier, result):

        txt = msg  # original schema validation error
        txt += f"\nRecord identifier {record_identifier} \nResult_identifier {result_identifier} \nReported result: {result}"
        super(SchemaValidationErrorDuringReport, self).__init__(txt)


class MissingConfigDataError(PipestatError):
    """Exception for invalid config file."""

    def __init__(self, msg):
        spacing = " " if msg[-1] in ["?", ".", "\n"] else "; "
        suggest = "For config format documentation please see: http://pipestat.databio.org/en/latest/db_config/"
        super(MissingConfigDataError, self).__init__(msg + spacing + suggest)


class PipestatStartupError(PipestatError):
    """Data error for local data associated with file backend"""

    def __init__(self, msg):
        super(PipestatStartupError, self).__init__(msg)


class PipestatDataError(PipestatError):
    """Data error for local data associated with file backend"""

    def __init__(self, msg):
        super(PipestatDataError, self).__init__(msg)


class PipestatSummarizeError(PipestatError):
    """Error for pipestat summarize"""

    def __init__(self, msg):
        super(PipestatSummarizeError, self).__init__(msg)


class InvalidTimeFormatError(PipestatError):
    """Data error for local data associated with file backend"""

    def __init__(self, msg):
        super(InvalidTimeFormatError, self).__init__(msg)


class PipestatDependencyError(PipestatError):
    """Dependency error"""

    def __init__(self, msg):
        super(PipestatDependencyError, self).__init__(msg)


class PipestatDatabaseError(PipestatError):
    """Database error"""

    def __init__(self, msg):
        super(PipestatDatabaseError, self).__init__(msg)


class PipestatPEPHubError(PipestatError):
    """PEPHub backend error"""

    def __init__(self, msg):
        super(PipestatPEPHubError, self).__init__(msg)


class InvalidTypeError(PipestatError):
    """Type of the reported value is not supported"""

    def __init__(self, type):
        super(InvalidTypeError, self).__init__(
            "'{}' is an invalid type. Only the following types are "
            "supported: {}".format(type, list(CLASSES_BY_TYPE.keys()))
        )


class IncompatibleClassError(PipestatError):
    """Class  of the reported value is not supported"""

    def __init__(self, cls, req_cls, type):
        super(IncompatibleClassError, self).__init__(
            "Incompatible value class for the declared result type ({}). "
            "Required: {}; got: {}".format(type, req_cls, cls)
        )


class UnrecognizedStatusError(PipestatError):
    """Exception for when a value to set as status isn't declared in the active status schema."""

    def __init__(self, status: str, known: Optional[Iterable[str]] = None):
        self._status = status
        msg = f"Unrecognized status: {status}"
        if known is not None:
            pass
        super(UnrecognizedStatusError, self).__init__(msg)

    @property
    def status(self):
        return self._status
