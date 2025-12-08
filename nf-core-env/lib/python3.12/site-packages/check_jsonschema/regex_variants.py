import enum
import re
import typing as t

import jsonschema
import regress


class RegexVariantName(enum.Enum):
    default = "default"
    nonunicode = "nonunicode"
    python = "python"


class RegexImplementation:
    """
    A high-level interface for getting at the different possible
    implementations of regex behaviors.
    """

    _concrete: "_ConcreteImplementation"

    def __init__(self, variant: RegexVariantName) -> None:
        self.variant = variant

        if self.variant == RegexVariantName.default:
            self._concrete = _RegressImplementation()
        elif self.variant == RegexVariantName.nonunicode:
            self._concrete = _NonunicodeRegressImplementation()
        else:
            self._concrete = _PythonImplementation()

        self.check_format = self._concrete.check_format
        self.pattern_keyword = self._concrete.pattern_keyword
        self.patternProperties_keyword = self._concrete.patternProperties_keyword


class _ConcreteImplementation(t.Protocol):
    def check_format(self, instance: t.Any) -> bool: ...

    def pattern_keyword(
        self, validator: t.Any, pattern: str, instance: str, schema: t.Any
    ) -> t.Iterator[jsonschema.ValidationError]: ...

    def patternProperties_keyword(
        self,
        validator: t.Any,
        patternProperties: dict[str, t.Any],
        instance: dict[str, t.Any],
        schema: t.Any,
    ) -> t.Iterator[jsonschema.ValidationError]: ...


class _RegressImplementation:
    def _compile_pattern(self, pattern: str) -> regress.Regex:
        return regress.Regex(pattern, flags="u")

    def check_format(self, instance: t.Any) -> bool:
        if not isinstance(instance, str):
            return True
        try:
            self._compile_pattern(instance)
        except regress.RegressError:
            return False
        return True

    def pattern_keyword(
        self, validator: t.Any, pattern: str, instance: str, schema: t.Any
    ) -> t.Iterator[jsonschema.ValidationError]:
        if not validator.is_type(instance, "string"):
            return

        regress_pattern = self._compile_pattern(pattern)
        if not regress_pattern.find(instance):
            yield jsonschema.ValidationError(f"{instance!r} does not match {pattern!r}")

    def patternProperties_keyword(
        self,
        validator: t.Any,
        patternProperties: dict[str, t.Any],
        instance: dict[str, t.Any],
        schema: t.Any,
    ) -> t.Iterator[jsonschema.ValidationError]:
        if not validator.is_type(instance, "object"):
            return

        for pattern, subschema in patternProperties.items():
            regress_pattern = self._compile_pattern(pattern)
            for k, v in instance.items():
                if regress_pattern.find(k):
                    yield from validator.descend(
                        v,
                        subschema,
                        path=k,
                        schema_path=pattern,
                    )


class _NonunicodeRegressImplementation(_RegressImplementation):
    def _compile_pattern(self, pattern: str) -> regress.Regex:
        return regress.Regex(pattern)


class _PythonImplementation:
    def check_format(self, instance: t.Any) -> bool:
        if not isinstance(instance, str):
            return True
        try:
            re.compile(instance)
        except re.error:
            return False
        return True

    def pattern_keyword(
        self, validator: t.Any, pattern: str, instance: str, schema: t.Any
    ) -> t.Iterator[jsonschema.ValidationError]:
        if not validator.is_type(instance, "string"):
            return

        re_pattern = re.compile(pattern)
        if not re_pattern.search(instance):
            yield jsonschema.ValidationError(f"{instance!r} does not match {pattern!r}")

    def patternProperties_keyword(
        self,
        validator: t.Any,
        patternProperties: dict[str, t.Any],
        instance: dict[str, t.Any],
        schema: t.Any,
    ) -> t.Iterator[jsonschema.ValidationError]:
        if not validator.is_type(instance, "object"):
            return

        for pattern, subschema in patternProperties.items():
            for k, v in instance.items():
                if re.search(pattern, k):
                    yield from validator.descend(
                        v,
                        subschema,
                        path=k,
                        schema_path=pattern,
                    )
