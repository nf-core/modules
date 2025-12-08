from __future__ import annotations

import copy

import jsonschema
import jsonschema.validators

from ..regex_variants import RegexImplementation
from .implementations import validate_rfc3339, validate_time

# all known format strings except for a selection from draft3 which have either
# been renamed or removed:
# - color
# - host-name
# - ip-address
KNOWN_FORMATS: tuple[str, ...] = (
    "date",
    "date-time",
    "duration",
    "email",
    "hostname",
    "idn-email",
    "idn-hostname",
    "ipv4",
    "ipv6",
    "iri",
    "iri-reference",
    "json-pointer",
    "regex",
    "relative-json-pointer",
    "time",
    "uri",
    "uri-reference",
    "uri-template",
    "uuid",
)


class FormatOptions:
    def __init__(
        self,
        *,
        regex_impl: RegexImplementation,
        enabled: bool = True,
        disabled_formats: tuple[str, ...] = (),
    ) -> None:
        self.enabled = enabled
        self.regex_impl = regex_impl
        self.disabled_formats = disabled_formats


def get_base_format_checker(schema_dialect: str | None) -> jsonschema.FormatChecker:
    # mypy does not consider a class whose instances match a protocol to match
    # `type[$PROTOCOL]` so this is considered a mismatch
    default_validator_cls: type[jsonschema.Validator] = (
        jsonschema.Draft202012Validator  # type:ignore[assignment]
    )
    # resolve the dialect, if given, to a validator class
    # default to the latest draft
    validator_class = jsonschema.validators.validator_for(
        {} if schema_dialect is None else {"$schema": schema_dialect},
        default=default_validator_cls,
    )
    return validator_class.FORMAT_CHECKER


def make_format_checker(
    opts: FormatOptions,
    schema_dialect: str | None = None,
) -> jsonschema.FormatChecker | None:
    if not opts.enabled:
        return None

    # customize around regex checking first
    checker = format_checker_for_regex_impl(opts.regex_impl)

    # add other custom format checks
    checker.checks("date-time")(validate_rfc3339)
    checker.checks("time")(validate_time)

    # remove the disabled checks, which may include the regex check
    for checkname in opts.disabled_formats:
        if checkname not in checker.checkers:
            continue
        del checker.checkers[checkname]

    return checker


def format_checker_for_regex_impl(
    regex_impl: RegexImplementation, schema_dialect: str | None = None
) -> jsonschema.FormatChecker:
    # convert to a schema-derived format checker, and copy it
    # for safe modification
    base_checker = get_base_format_checker(schema_dialect)
    checker = copy.deepcopy(base_checker)

    # replace the regex check
    del checker.checkers["regex"]
    checker.checks("regex")(regex_impl.check_format)

    return checker
