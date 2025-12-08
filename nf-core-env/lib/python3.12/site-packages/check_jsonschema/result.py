from __future__ import annotations

import pathlib

import jsonschema

from .parsers import ParseError


class CheckResult:
    def __init__(self) -> None:
        self.validation_errors: dict[str, list[jsonschema.ValidationError]] = {}
        self.parse_errors: dict[str, list[ParseError]] = {}
        self.successes: list[str] = []

    @property
    def success(self) -> bool:
        return not (bool(self.parse_errors) or bool(self.validation_errors))

    def record_validation_success(self, path: pathlib.Path | str) -> None:
        self.successes.append(str(path))

    def record_validation_error(
        self, path: pathlib.Path | str, err: jsonschema.ValidationError
    ) -> None:
        filename = str(path)
        if filename not in self.validation_errors:
            self.validation_errors[filename] = []
        self.validation_errors[filename].append(err)

    def record_parse_error(self, path: pathlib.Path | str, err: ParseError) -> None:
        filename = str(path)
        if filename not in self.parse_errors:
            self.parse_errors[filename] = []
        self.parse_errors[filename].append(err)
