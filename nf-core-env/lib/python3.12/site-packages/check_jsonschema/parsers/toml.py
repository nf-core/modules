from __future__ import annotations

import datetime
import sys
import typing as t

if sys.version_info < (3, 11):
    import tomli as toml_implementation
else:
    import tomllib as toml_implementation


def _normalize(data: t.Any) -> t.Any:
    """
    Normalize TOML data to fit the requirements to be JSON-encodeable.

    Currently this applies the following transformations:

        offset-aware datetime.datetime values are converted to strings using isoformat()
        naive datetime.datetime values are converted to strings using isoformat() + "Z"

        offset-aware datetime.time values are converted to strings using isoformat()
        naive datetime.time values are converted to strings using isoformat() + "Z"

        datetime.date values are converted to strings using isoformat()
    """
    if isinstance(data, dict):
        return {k: _normalize(v) for k, v in data.items()}
    elif isinstance(data, list):
        return [_normalize(x) for x in data]
    else:
        # python's datetime will format to an ISO partial time when handling a naive
        # time/datetime , but JSON Schema format validation specifies that date-time is
        # taken from RFC3339, which defines "date-time" as including 'Z|offset'
        # the specification for "time" is less clear because JSON Schema does not specify
        # which RFC3339 definition should be used, and the RFC has no format named "time",
        # only "full-time" (with Z|offset) and "partial-time" (no offset)
        #
        # rfc3339_validator (used by 'jsonschema') requires the offset, so we will do the
        # same
        if isinstance(data, datetime.datetime) or isinstance(data, datetime.time):
            if data.tzinfo is None:
                return data.isoformat() + "Z"
            return data.isoformat()
        elif isinstance(data, datetime.date):
            return data.isoformat()
        return data


ParseError: type[Exception] = toml_implementation.TOMLDecodeError


def load(stream: t.IO[bytes]) -> t.Any:
    data = toml_implementation.load(stream)
    return _normalize(data)
