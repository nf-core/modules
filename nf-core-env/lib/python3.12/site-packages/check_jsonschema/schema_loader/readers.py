from __future__ import annotations

import io
import json
import sys
import typing as t

import ruamel.yaml

from ..cachedownloader import CacheDownloader
from ..parsers import ParseError, ParserSet
from ..utils import filename2path
from .errors import SchemaParseError

yaml = ruamel.yaml.YAML(typ="safe")


class _UnsetType:
    pass


_UNSET = _UnsetType()


def _run_load_callback(schema_location: str, callback: t.Callable) -> dict:
    try:
        schema = callback()
    # only local loads can raise the YAMLError, but catch for both cases for simplicity
    except (ValueError, ruamel.yaml.error.YAMLError) as e:
        raise SchemaParseError(schema_location) from e
    if not isinstance(schema, dict):
        raise SchemaParseError(schema_location)
    return schema


class LocalSchemaReader:
    def __init__(self, filename: str) -> None:
        self.path = filename2path(filename)
        self.filename = str(self.path)
        self.parsers = ParserSet()
        self._parsed_schema: dict | _UnsetType = _UNSET

    def get_retrieval_uri(self) -> str | None:
        return self.path.as_uri()

    def _read_impl(self) -> t.Any:
        return self.parsers.parse_file(self.path, default_filetype="json")

    def read_schema(self) -> dict:
        if self._parsed_schema is _UNSET:
            self._parsed_schema = _run_load_callback(self.filename, self._read_impl)
        return t.cast(dict, self._parsed_schema)


class StdinSchemaReader:
    def __init__(self) -> None:
        self.parsers = ParserSet()
        self._parsed_schema: dict | _UnsetType = _UNSET

    def get_retrieval_uri(self) -> str | None:
        return None

    def read_schema(self) -> dict:
        if self._parsed_schema is _UNSET:
            try:
                self._parsed_schema = json.load(sys.stdin)
            except ValueError as e:
                raise ParseError("Failed to parse JSON from stdin") from e
        return t.cast(dict, self._parsed_schema)


class HttpSchemaReader:
    def __init__(
        self,
        url: str,
        disable_cache: bool,
    ) -> None:
        self.url = url
        self.parsers = ParserSet()
        self.downloader = CacheDownloader("schemas", disable_cache=disable_cache).bind(
            url, validation_callback=self._parse
        )
        self._parsed_schema: dict | _UnsetType = _UNSET

    def _parse(self, schema_bytes: bytes) -> t.Any:
        return self.parsers.parse_data_with_path(
            io.BytesIO(schema_bytes), self.url, default_filetype="json"
        )

    def get_retrieval_uri(self) -> str | None:
        return self.url

    def _read_impl(self) -> t.Any:
        with self.downloader.open() as fp:
            return self._parse(fp.read())

    def read_schema(self) -> dict:
        if self._parsed_schema is _UNSET:
            self._parsed_schema = _run_load_callback(self.url, self._read_impl)
        return t.cast(dict, self._parsed_schema)
