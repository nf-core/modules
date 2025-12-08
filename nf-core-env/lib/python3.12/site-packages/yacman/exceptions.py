"""Package exception types"""

__all__ = ["FileFormatError", "AliasError", "UndefinedAliasError"]


class FileFormatError(Exception):
    """Exception for invalid file format."""

    pass


class AliasError(Exception):
    """Alias related error."""

    pass


class UndefinedAliasError(AliasError):
    """Alias is is not defined."""

    pass
