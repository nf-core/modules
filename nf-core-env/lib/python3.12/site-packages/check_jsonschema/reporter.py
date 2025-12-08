"""
Output formatters are called "reporters" because they take the result of validation
and report it back to the user.
"""

from __future__ import annotations

import abc
import json
import textwrap
import typing as t

import click
import jsonschema

from . import format_errors
from .parsers import ParseError
from .result import CheckResult
from .utils import iter_validation_error


class Reporter(abc.ABC):
    def __init__(self, *, verbosity: int, **kwargs: t.Any) -> None:
        self.verbosity = verbosity
        super().__init__(**kwargs)

    @abc.abstractmethod
    def report_success(self, result: CheckResult) -> None:
        raise NotImplementedError

    @abc.abstractmethod
    def report_errors(self, result: CheckResult) -> None:
        raise NotImplementedError

    def report_result(self, result: CheckResult) -> None:
        if result.success:
            self.report_success(result)
        else:
            self.report_errors(result)


class TextReporter(Reporter):
    def __init__(
        self,
        *,
        verbosity: int,
        stream: t.TextIO | None = None,  # default stream is stdout (None)
    ) -> None:
        super().__init__(verbosity=verbosity)
        self.stream = stream

    def _echo(self, s: str, *, indent: int = 0) -> None:
        click.echo(" " * indent + s, file=self.stream)

    def report_success(self, result: CheckResult) -> None:
        if self.verbosity < 1:
            return
        ok = click.style("ok", fg="green")
        self._echo(f"{ok} -- validation done")
        if self.verbosity > 1:
            self._echo("The following files were checked:")
            for filename in result.successes:
                self._echo(f"  {filename}")

    def _format_validation_error_message(
        self, err: jsonschema.ValidationError, filename: str | None = None
    ) -> str:
        error_loc = err.json_path
        if filename:
            error_loc = f"{filename}::{error_loc}"
        error_loc = click.style(error_loc, fg="yellow")
        return f"{error_loc}: {err.message}"

    def _show_validation_error(
        self,
        filename: str,
        err: jsonschema.ValidationError,
    ) -> None:
        self._echo(
            self._format_validation_error_message(err, filename=filename), indent=2
        )
        if err.context:
            best_match = jsonschema.exceptions.best_match(err.context)
            self._echo("Underlying errors caused this.", indent=2)
            self._echo("")
            self._echo("Best Match:", indent=2)
            self._echo(self._format_validation_error_message(best_match), indent=4)

            best_deep_match = find_best_deep_match(err)
            if best_deep_match != best_match:
                self._echo("Best Deep Match:", indent=2)
                self._echo(
                    self._format_validation_error_message(best_deep_match), indent=4
                )

            if self.verbosity > 1:
                self._echo("All Errors:", indent=2)
                for e in iter_validation_error(err):
                    self._echo(self._format_validation_error_message(e), indent=4)
            else:
                num_other_errors = len(list(iter_validation_error(err))) - 1
                if best_deep_match != best_match:
                    num_other_errors -= 1
                if num_other_errors > 0:
                    self._echo("")
                    self._echo(
                        f"{click.style(str(num_other_errors), fg='yellow')} other "
                        "errors were produced. "
                        "Use '--verbose' to see all errors.",
                        indent=2,
                    )

    def _show_parse_error(self, filename: str, err: ParseError) -> None:
        if self.verbosity < 2:
            mode: t.Literal["minimal", "short", "full"] = "minimal"
        elif self.verbosity < 3:
            mode = "short"
        else:
            mode = "full"
        self._echo(textwrap.indent(format_errors.format_error(err, mode=mode), "  "))

    def report_errors(self, result: CheckResult) -> None:
        if self.verbosity < 1:
            return
        if result.parse_errors:
            self._echo("Several files failed to parse.")
            for filename, errors in result.parse_errors.items():
                for err in errors:
                    self._show_parse_error(filename, err)
        if result.validation_errors:
            self._echo("Schema validation errors were encountered.")
            for filename, parse_errors in result.validation_errors.items():
                for parse_err in parse_errors:
                    self._show_validation_error(filename, parse_err)


class JsonReporter(Reporter):
    def __init__(self, *, verbosity: int, pretty: bool = True) -> None:
        super().__init__(verbosity=verbosity)
        # default to pretty output, can add a switch to disable this in the future
        self.pretty = pretty

    def _dump(self, data: t.Any) -> None:
        if self.pretty:
            click.echo(json.dumps(data, indent=2, separators=(",", ": ")))
        else:
            click.echo(json.dumps(data, separators=(",", ":")))

    def report_success(self, result: CheckResult) -> None:
        report_obj: dict[str, t.Any] = {"status": "ok"}
        if self.verbosity > 0:
            report_obj["errors"] = []
        if self.verbosity > 1:
            report_obj["checked_paths"] = list(result.successes)
        self._dump(report_obj)

    def _dump_error_map(
        self,
        error_map: dict[str, list[jsonschema.ValidationError]],
    ) -> t.Iterator[dict]:
        for filename, errors in error_map.items():
            for err in errors:
                item = {
                    "filename": filename,
                    "path": err.json_path,
                    "message": err.message,
                    "has_sub_errors": bool(err.context),
                }
                if err.context:
                    best_match = jsonschema.exceptions.best_match(err.context)
                    best_deep_match = find_best_deep_match(err)
                    item["best_match"] = {
                        "path": best_match.json_path,
                        "message": best_match.message,
                    }
                    item["best_deep_match"] = {
                        "path": best_deep_match.json_path,
                        "message": best_deep_match.message,
                    }
                    num_sub_errors = len(list(iter_validation_error(err))) - 1
                    item["num_sub_errors"] = num_sub_errors
                    if self.verbosity > 1:
                        item["sub_errors"] = [
                            {"path": suberr.json_path, "message": suberr.message}
                            for suberr in iter_validation_error(err)
                        ]

                yield item

    def _dump_parse_errors(
        self,
        error_map: dict[str, list[ParseError]],
    ) -> t.Iterator[dict]:
        for filename, errors in error_map.items():
            for err in errors:
                yield {
                    "filename": filename,
                    "message": str(err),
                }

    def report_errors(self, result: CheckResult) -> None:
        report_obj: dict[str, t.Any] = {"status": "fail"}
        if self.verbosity > 1:
            report_obj["successes"] = list(result.successes)
        if self.verbosity > 0:
            report_obj["errors"] = list(self._dump_error_map(result.validation_errors))
            report_obj["parse_errors"] = list(
                self._dump_parse_errors(result.parse_errors)
            )
        self._dump(report_obj)


REPORTER_BY_NAME: dict[str, type[Reporter]] = {
    "text": TextReporter,
    "json": JsonReporter,
}


def _deep_match_relevance(error: jsonschema.ValidationError) -> tuple[bool | int, ...]:
    validator = error.validator
    return (
        validator not in ("anyOf", "oneOf"),
        len(error.absolute_path),
        -len(error.path),
    )


def find_best_deep_match(
    errors: jsonschema.ValidationError,
) -> jsonschema.ValidationError:
    return max(iter_validation_error(errors), key=_deep_match_relevance)
