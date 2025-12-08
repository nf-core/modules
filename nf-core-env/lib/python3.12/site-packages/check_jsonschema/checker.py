from __future__ import annotations

import pathlib
import typing as t

import click
import jsonschema
import referencing.exceptions

from . import format_errors
from .formats import FormatOptions
from .instance_loader import InstanceLoader
from .parsers import ParseError
from .regex_variants import RegexImplementation
from .reporter import Reporter
from .result import CheckResult
from .schema_loader import SchemaLoaderBase, SchemaParseError, UnsupportedUrlScheme


class _Exit(Exception):
    def __init__(self, code: int) -> None:
        self.code = code


class SchemaChecker:
    def __init__(
        self,
        schema_loader: SchemaLoaderBase,
        instance_loader: InstanceLoader,
        reporter: Reporter,
        *,
        format_opts: FormatOptions,
        regex_impl: RegexImplementation,
        traceback_mode: t.Literal["minimal", "short", "full"] = "short",
        fill_defaults: bool = False,
    ) -> None:
        self._schema_loader = schema_loader
        self._instance_loader = instance_loader
        self._reporter = reporter

        self._format_opts = format_opts
        self._regex_impl = regex_impl
        self._traceback_mode = traceback_mode
        self._fill_defaults = fill_defaults

    def _fail(self, msg: str, err: Exception | None = None) -> t.NoReturn:
        click.echo(msg, err=True)
        if err is not None:
            format_errors.print_error(err, mode=self._traceback_mode)
        raise _Exit(1)

    def get_validator(
        self, path: pathlib.Path | str, doc: dict[str, t.Any]
    ) -> jsonschema.protocols.Validator:
        try:
            return self._schema_loader.get_validator(
                path, doc, self._format_opts, self._regex_impl, self._fill_defaults
            )
        except SchemaParseError as e:
            self._fail("Error: schemafile could not be parsed as JSON", e)
        except jsonschema.SchemaError as e:
            self._fail("Error: schemafile was not valid\n", e)
        except UnsupportedUrlScheme as e:
            self._fail(f"Error: {e}\n", e)
        except Exception as e:
            self._fail("Error: Unexpected Error building schema validator", e)

    def _build_result(self) -> CheckResult:
        result = CheckResult()
        for path, data in self._instance_loader.iter_files():
            if isinstance(data, ParseError):
                result.record_parse_error(path, data)
            else:
                validator = self.get_validator(path, data)
                passing = True
                for err in validator.iter_errors(data):
                    result.record_validation_error(path, err)
                    passing = False
                if passing:
                    result.record_validation_success(path)
        return result

    def _run(self) -> None:
        try:
            result = self._build_result()
        except (
            referencing.exceptions.NoSuchResource,
            referencing.exceptions.Unretrievable,
            referencing.exceptions.Unresolvable,
        ) as e:
            self._fail("Failure resolving $ref within schema\n", e)

        self._reporter.report_result(result)
        if not result.success:
            raise _Exit(1)

    def run(self) -> int:
        try:
            self._run()
        except _Exit as e:
            return e.code
        return 0
