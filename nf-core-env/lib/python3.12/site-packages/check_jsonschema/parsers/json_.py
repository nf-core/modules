from __future__ import annotations

import json
import typing as t

try:
    import orjson

    has_orjson = True
except ImportError:
    has_orjson = False

JSONDecodeError = json.JSONDecodeError


def load(stream: t.IO[bytes]) -> t.Any:
    bin_data = stream.read()
    # if orjson is available, try it first
    if has_orjson:
        # in the event of a decode error, it may be that the data contains
        # `Infinity`, `-Infinity`, or `NaN`
        #
        # however, do not failover to stdlib JSON -- it is not obvious that there's any
        # need for check-jsonschema to support these invalid JSON datatypes
        # if users encounter issues with this behavior in the future, we can revisit how
        # JSON loading is handled
        return orjson.loads(bin_data)
    # failover to stdlib json
    return json.loads(bin_data)
