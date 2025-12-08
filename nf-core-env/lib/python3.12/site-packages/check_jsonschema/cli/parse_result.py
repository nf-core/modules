from __future__ import annotations

import enum
import typing as t

import click
import jsonschema

from ..formats import FormatOptions
from ..regex_variants import RegexImplementation, RegexVariantName
from ..transforms import Transform


class SchemaLoadingMode(enum.Enum):
    filepath = "filepath"
    builtin = "builtin"
    metaschema = "metaschema"


class ParseResult:
    def __init__(self) -> None:
        # primary options: schema + instances
        self.schema_mode: SchemaLoadingMode = SchemaLoadingMode.filepath
        self.schema_path: str | None = None
        self.base_uri: str | None = None
        self.instancefiles: tuple[t.IO[bytes], ...] = ()
        # cache controls
        self.disable_cache: bool = False
        self.cache_filename: str | None = None
        # filetype detection (JSON, YAML, TOML, etc)
        self.default_filetype: str = "json"
        self.force_filetype: str | None = None
        # data-transform (for Azure Pipelines and potentially future transforms)
        self.data_transform: Transform | None = None
        # validation behavioral controls
        self.validator_class: type[jsonschema.protocols.Validator] | None = None
        self.fill_defaults: bool = False
        # regex format options
        self.disable_all_formats: bool = False
        self.disable_formats: tuple[str, ...] = ()
        self.regex_variant: RegexVariantName = RegexVariantName.default
        # error and output controls
        self.verbosity: int = 1
        self.traceback_mode: t.Literal["short", "full"] = "short"
        self.output_format: str = "text"

    def set_regex_variant(
        self,
        variant_opt: t.Literal["python", "nonunicode", "default"] | None,
        *,
        legacy_opt: t.Literal["python", "nonunicode", "default"] | None = None,
    ) -> None:
        variant_name: t.Literal["python", "nonunicode", "default"] | None = (
            variant_opt or legacy_opt
        )
        if variant_name:
            self.regex_variant = RegexVariantName(variant_name)

    def set_schema(
        self, schemafile: str | None, builtin_schema: str | None, check_metaschema: bool
    ) -> None:
        mutex_arg_count = sum(
            1 if x else 0 for x in (schemafile, builtin_schema, check_metaschema)
        )
        if mutex_arg_count == 0:
            raise click.UsageError(
                "Either --schemafile, --builtin-schema, or --check-metaschema "
                "must be provided"
            )
        if mutex_arg_count > 1:
            raise click.UsageError(
                "--schemafile, --builtin-schema, and --check-metaschema "
                "are mutually exclusive"
            )

        if schemafile:
            self.schema_mode = SchemaLoadingMode.filepath
            self.schema_path = schemafile
        elif builtin_schema:
            self.schema_mode = SchemaLoadingMode.builtin
            self.schema_path = builtin_schema
        else:
            self.schema_mode = SchemaLoadingMode.metaschema

    def set_validator(
        self, validator_class: type[jsonschema.protocols.Validator] | None
    ) -> None:
        if validator_class is None:
            return
        if self.schema_mode != SchemaLoadingMode.filepath:
            raise click.UsageError(
                "--validator-class can only be used with --schemafile for schema loading"
            )
        self.validator_class = validator_class

    @property
    def format_opts(self) -> FormatOptions:
        return FormatOptions(
            regex_impl=RegexImplementation(self.regex_variant),
            enabled=not self.disable_all_formats,
            disabled_formats=self.disable_formats,
        )
