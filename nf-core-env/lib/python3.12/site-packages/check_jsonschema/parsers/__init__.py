from __future__ import annotations

import io
import pathlib
import typing as t

import ruamel.yaml

from ..identify_filetype import path_to_type
from . import json5, json_, toml, yaml

_PARSER_ERRORS: set[type[Exception]] = {
    json_.JSONDecodeError,
    yaml.ParseError,
    toml.ParseError,
}
DEFAULT_LOAD_FUNC_BY_TAG: dict[str, t.Callable[[t.IO[bytes]], t.Any]] = {
    "json": json_.load,
    "toml": toml.load,
}
SUPPORTED_FILE_FORMATS = ["json", "toml", "yaml"]
if json5.ENABLED:
    SUPPORTED_FILE_FORMATS.append("json5")
    DEFAULT_LOAD_FUNC_BY_TAG["json5"] = json5.load
    _PARSER_ERRORS.add(json5.ParseError)
MISSING_SUPPORT_MESSAGES: dict[str, str] = {
    "json5": json5.MISSING_SUPPORT_MESSAGE,
}
LOADING_FAILURE_ERROR_TYPES: tuple[type[Exception], ...] = tuple(_PARSER_ERRORS)


class ParseError(ValueError):
    pass


class BadFileTypeError(ParseError):
    pass


class FailedFileLoadError(ParseError):
    pass


class ParserSet:
    def __init__(
        self,
        *,
        modify_yaml_implementation: t.Callable[[ruamel.yaml.YAML], None] | None = None,
        supported_formats: t.Sequence[str] | None = None,
    ) -> None:
        yaml_impl = yaml.construct_yaml_implementation()
        failover_yaml_impl = yaml.construct_yaml_implementation(pure=True)
        if modify_yaml_implementation:
            modify_yaml_implementation(yaml_impl)
            modify_yaml_implementation(failover_yaml_impl)
        base_by_tag = {
            "yaml": yaml.impl2loader(yaml_impl, failover_yaml_impl),
            **DEFAULT_LOAD_FUNC_BY_TAG,
        }
        if supported_formats is None:
            self._by_tag = base_by_tag
        else:
            self._by_tag = {
                k: v for k, v in base_by_tag.items() if k in supported_formats
            }

    def get(
        self,
        path: pathlib.Path | str,
        default_filetype: str,
        force_filetype: str | None = None,
    ) -> t.Callable[[t.IO[bytes]], t.Any]:
        if force_filetype:
            filetype = force_filetype
        else:
            filetype = path_to_type(path, default_type=default_filetype)

        if filetype in self._by_tag:
            return self._by_tag[filetype]

        if filetype in MISSING_SUPPORT_MESSAGES:
            raise BadFileTypeError(
                f"cannot parse {path} because support is missing for {filetype}\n"
                + MISSING_SUPPORT_MESSAGES[filetype]
            )
        raise BadFileTypeError(
            f"cannot parse {path} as it is not one of the supported filetypes: "
            + ",".join(self._by_tag.keys())
        )

    def parse_data_with_path(
        self,
        data: t.IO[bytes] | bytes,
        path: pathlib.Path | str,
        default_filetype: str,
        force_filetype: str | None = None,
    ) -> t.Any:
        loadfunc = self.get(path, default_filetype, force_filetype)
        try:
            if isinstance(data, bytes):
                data = io.BytesIO(data)
            return loadfunc(data)
        except LOADING_FAILURE_ERROR_TYPES as e:
            raise FailedFileLoadError(f"Failed to parse {path}") from e

    def parse_file(
        self,
        path: pathlib.Path | str,
        default_filetype: str,
        force_filetype: str | None = None,
    ) -> t.Any:
        with open(path, "rb") as fp:
            return self.parse_data_with_path(fp, path, default_filetype, force_filetype)
