from __future__ import annotations

import os
import textwrap
import typing as t

import click
import jsonschema

from ..catalog import CUSTOM_SCHEMA_NAMES, SCHEMA_CATALOG
from ..checker import SchemaChecker
from ..formats import KNOWN_FORMATS
from ..instance_loader import InstanceLoader
from ..parsers import SUPPORTED_FILE_FORMATS
from ..regex_variants import RegexImplementation, RegexVariantName
from ..reporter import REPORTER_BY_NAME, Reporter
from ..schema_loader import (
    BuiltinSchemaLoader,
    MetaSchemaLoader,
    SchemaLoader,
    SchemaLoaderBase,
)
from ..transforms import TRANSFORM_LIBRARY
from .param_types import CommaDelimitedList, LazyBinaryReadFile, ValidatorClassName
from .parse_result import ParseResult, SchemaLoadingMode

BUILTIN_SCHEMA_NAMES = [f"vendor.{k}" for k in SCHEMA_CATALOG.keys()] + [
    f"custom.{k}" for k in CUSTOM_SCHEMA_NAMES
]
BUILTIN_SCHEMA_CHOICES = (
    BUILTIN_SCHEMA_NAMES + list(SCHEMA_CATALOG.keys()) + CUSTOM_SCHEMA_NAMES
)


def set_color_mode(ctx: click.Context, param: str, value: str) -> None:
    if "NO_COLOR" in os.environ:
        ctx.color = False
    else:
        ctx.color = {
            "auto": None,
            "always": True,
            "never": False,
        }[value]


def pretty_helptext_list(values: list[str] | tuple[str, ...]) -> str:
    return textwrap.indent(
        "\n".join(
            textwrap.wrap(
                ", ".join(values),
                width=75,
                break_long_words=False,
                break_on_hyphens=False,
            ),
        ),
        "    ",
    )


@click.command(
    "check-jsonschema",
    help="""\
Check JSON and YAML files against a JSON Schema.

The schema is specified either with '--schemafile' or with '--builtin-schema'.

'check-jsonschema' supports format checks with appropriate libraries installed,
including the following formats by default:
    date, date-time, email, ipv4, ipv6, regex, uuid

\b
For handling of regexes, there are multiple modes which can be specified with
'--regex-variant':
    default    |  use ECMAScript regex syntax (via regress)
    nonunicode |  use ECMAScript regex syntax, but in non-unicode mode (via regress)
    python     |  use python regex syntax

\b
The '--builtin-schema' flag supports the following schema names:
"""
    + pretty_helptext_list(BUILTIN_SCHEMA_NAMES)
    + """\

\b
The '--disable-formats' flag supports the following formats:
"""
    + pretty_helptext_list(KNOWN_FORMATS),
)
@click.help_option("-h", "--help")
@click.version_option()
@click.option(
    "--schemafile",
    help=(
        "The path to a file containing the JSON Schema to use or an "
        "HTTP(S) URI for the schema. If a remote file is used, "
        "it will be downloaded and cached locally based on mtime. "
        "Use '-' for stdin."
    ),
    metavar="[PATH|URI]",
)
@click.option(
    "--base-uri",
    help=(
        "Override the base URI for the schema. The default behavior is to "
        "follow the behavior specified by the JSON Schema spec, which is to "
        "prefer an explicit '$id' and failover to the retrieval URI."
    ),
)
@click.option(
    "--builtin-schema",
    help="The name of an internal schema to use for '--schemafile'",
    type=click.Choice(BUILTIN_SCHEMA_CHOICES, case_sensitive=False),
    metavar="BUILTIN_SCHEMA_NAME",
)
@click.option(
    "--check-metaschema",
    is_flag=True,
    help=(
        "Instead of validating the instances against a schema, treat each file as a "
        "schema and validate them under their matching metaschemas."
    ),
)
@click.option(
    "--no-cache",
    is_flag=True,
    help="Disable schema caching. Always download remote schemas.",
)
@click.option(
    "--cache-filename", help="Deprecated. This option no longer has any effect."
)
@click.option(
    "--disable-formats",
    multiple=True,
    help=(
        "Disable specific format checks in the schema. "
        "Pass '*' to disable all format checks."
    ),
    type=CommaDelimitedList(choices=("*", *KNOWN_FORMATS)),
    metavar="{*|FORMAT,FORMAT,...}",
)
@click.option(
    "--format-regex",
    hidden=True,
    help="Legacy name for `--regex-variant`.",
    default=None,
    type=click.Choice([x.value for x in RegexVariantName], case_sensitive=False),
)
@click.option(
    "--regex-variant",
    help=(
        "Name of which regex dialect should be used for format checking "
        "and 'pattern' matching."
    ),
    default=None,
    type=click.Choice([x.value for x in RegexVariantName], case_sensitive=False),
)
@click.option(
    "--default-filetype",
    help="A default filetype to assume when a file's type is not detected",
    default="json",
    show_default=True,
    type=click.Choice(SUPPORTED_FILE_FORMATS, case_sensitive=True),
)
@click.option(
    "--force-filetype",
    help="Force a file type to use when parsing instance files",
    type=click.Choice(SUPPORTED_FILE_FORMATS, case_sensitive=True),
)
@click.option(
    "--traceback-mode",
    help=(
        "Set the mode of presentation for error traces. "
        "Defaults to shortened tracebacks."
    ),
    type=click.Choice(("full", "short")),
    default="short",
)
@click.option(
    "--data-transform",
    help=(
        "Select a builtin transform which should be applied to instancefiles before "
        "they are checked."
    ),
    type=click.Choice(tuple(TRANSFORM_LIBRARY.keys())),
)
@click.option(
    "--fill-defaults",
    help=(
        "Autofill 'default' values prior to validation. "
        "This may conflict with certain third-party validators used with "
        "'--validator-class'"
    ),
    is_flag=True,
)
@click.option(
    "--validator-class",
    help=(
        "The fully qualified name of a python validator to use in place of "
        "the 'jsonschema' library validators, in the form of '<package>:<class>'. "
        "The validator must be importable in the same environment where "
        "'check-jsonschema' is run."
    ),
    type=ValidatorClassName(),
)
@click.option(
    "-o",
    "--output-format",
    help="Which output format to use.",
    type=click.Choice(tuple(REPORTER_BY_NAME.keys()), case_sensitive=False),
    default="text",
)
@click.option(
    "--color",
    help="Force or disable colorized output. Defaults to autodetection.",
    default="auto",
    type=click.Choice(("auto", "always", "never")),
    callback=set_color_mode,
    expose_value=False,
)
@click.option(
    "-v",
    "--verbose",
    help=(
        "Increase output verbosity. On validation errors, this may be especially "
        "useful when oneOf or anyOf is used in the schema."
    ),
    count=True,
)
@click.option(
    "-q",
    "--quiet",
    help="Reduce output verbosity",
    count=True,
)
@click.argument(
    "instancefiles", required=True, nargs=-1, type=LazyBinaryReadFile("rb", lazy=True)
)
def main(
    *,
    schemafile: str | None,
    builtin_schema: str | None,
    base_uri: str | None,
    check_metaschema: bool,
    no_cache: bool,
    cache_filename: str | None,
    disable_formats: tuple[list[str], ...],
    format_regex: t.Literal["python", "nonunicode", "default"] | None,
    regex_variant: t.Literal["python", "nonunicode", "default"] | None,
    default_filetype: t.Literal["json", "yaml", "toml", "json5"],
    force_filetype: t.Literal["json", "yaml", "toml", "json5"] | None,
    traceback_mode: t.Literal["full", "short"],
    data_transform: t.Literal["azure-pipelines", "gitlab-ci"] | None,
    fill_defaults: bool,
    validator_class: type[jsonschema.protocols.Validator] | None,
    output_format: t.Literal["text", "json"],
    verbose: int,
    quiet: int,
    instancefiles: tuple[t.IO[bytes], ...],
) -> None:
    args = ParseResult()

    args.set_regex_variant(regex_variant, legacy_opt=format_regex)

    args.set_schema(schemafile, builtin_schema, check_metaschema)
    args.set_validator(validator_class)

    args.base_uri = base_uri
    args.instancefiles = instancefiles

    normalized_disable_formats: tuple[str, ...] = tuple(
        f for sublist in disable_formats for f in sublist
    )
    if "*" in normalized_disable_formats:
        args.disable_all_formats = True
    else:
        args.disable_formats = normalized_disable_formats

    args.disable_cache = no_cache
    args.default_filetype = default_filetype
    args.force_filetype = force_filetype
    args.fill_defaults = fill_defaults
    if data_transform is not None:
        args.data_transform = TRANSFORM_LIBRARY[data_transform]

    # verbosity behavior:
    # - default is 1
    # - count '-v'
    # - subtract count of '-q'
    args.verbosity = 1 + verbose - quiet
    args.traceback_mode = traceback_mode
    args.output_format = output_format

    execute(args)


# separate parsing from execution for simpler mocking for unit tests


def build_schema_loader(args: ParseResult) -> SchemaLoaderBase:
    if args.schema_mode == SchemaLoadingMode.metaschema:
        return MetaSchemaLoader(base_uri=args.base_uri)
    elif args.schema_mode == SchemaLoadingMode.builtin:
        assert args.schema_path is not None
        return BuiltinSchemaLoader(args.schema_path, base_uri=args.base_uri)
    elif args.schema_mode == SchemaLoadingMode.filepath:
        assert args.schema_path is not None
        return SchemaLoader(
            args.schema_path,
            disable_cache=args.disable_cache,
            base_uri=args.base_uri,
            validator_class=args.validator_class,
        )
    else:
        raise NotImplementedError("no valid schema option provided")


def build_instance_loader(args: ParseResult) -> InstanceLoader:
    return InstanceLoader(
        args.instancefiles,
        default_filetype=args.default_filetype,
        force_filetype=args.force_filetype,
        data_transform=args.data_transform,
    )


def build_reporter(args: ParseResult) -> Reporter:
    cls = REPORTER_BY_NAME[args.output_format]
    return cls(verbosity=args.verbosity)


def build_checker(args: ParseResult) -> SchemaChecker:
    schema_loader = build_schema_loader(args)
    instance_loader = build_instance_loader(args)
    reporter = build_reporter(args)
    return SchemaChecker(
        schema_loader,
        instance_loader,
        reporter,
        format_opts=args.format_opts,
        regex_impl=RegexImplementation(args.regex_variant),
        traceback_mode=args.traceback_mode,
        fill_defaults=args.fill_defaults,
    )


def execute(args: ParseResult) -> None:
    checker = build_checker(args)
    ret = checker.run()
    click.get_current_context().exit(ret)
