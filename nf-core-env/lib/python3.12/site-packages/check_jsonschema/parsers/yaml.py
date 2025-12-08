from __future__ import annotations

import typing as t
import warnings

import ruamel.yaml

ParseError = ruamel.yaml.YAMLError


def construct_yaml_implementation(
    typ: str = "safe", pure: bool = False
) -> ruamel.yaml.YAML:
    implementation = ruamel.yaml.YAML(typ=typ, pure=pure)

    # workaround global state
    # see: https://sourceforge.net/p/ruamel-yaml/tickets/341/
    class GeneratedSafeConstructor(ruamel.yaml.SafeConstructor):
        pass

    implementation.Constructor = GeneratedSafeConstructor

    # ruamel.yaml parses timestamp values into datetime.datetime values
    # however, JSON does not support native datetimes, so JSON Schema formats for
    # dates apply to strings
    # Turn off this feature, instructing the parser to load datetimes as strings
    implementation.constructor.yaml_constructors["tag:yaml.org,2002:timestamp"] = (
        implementation.constructor.yaml_constructors["tag:yaml.org,2002:str"]
    )

    return implementation


def _normalize(data: t.Any) -> t.Any:
    """
    Normalize YAML data to fit the requirements to be JSON-encodeable.

    Currently this applies the following transformation:
        dict keys are converted to strings

    Additional tweaks can be added in this layer in the future if necessary.
    """
    if isinstance(data, dict):
        return {str(k): _normalize(v) for k, v in data.items()}
    elif isinstance(data, list):
        return [_normalize(x) for x in data]
    else:
        return data


_data_sentinel = object()


def impl2loader(
    primary: ruamel.yaml.YAML, *fallbacks: ruamel.yaml.YAML
) -> t.Callable[[t.IO[bytes]], t.Any]:
    def load(stream: t.IO[bytes]) -> t.Any:
        stream_bytes = stream.read()
        lasterr: ruamel.yaml.YAMLError | None = None
        data: t.Any = _data_sentinel
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ruamel.yaml.error.ReusedAnchorWarning)
            for impl in [primary] + list(fallbacks):
                try:
                    data = impl.load(stream_bytes)
                except ruamel.yaml.YAMLError as e:
                    lasterr = e
                else:
                    break
        if data is _data_sentinel and lasterr is not None:
            raise lasterr
        return _normalize(data)

    return load
