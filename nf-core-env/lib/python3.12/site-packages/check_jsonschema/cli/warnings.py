from __future__ import annotations

import typing as t
import warnings

import click


def deprecation_warning_callback(
    optstring: str, *, is_flag: bool = False, append_message: str | None = None
) -> t.Callable[[click.Context, click.Parameter, t.Any], t.Any]:
    def callback(ctx: click.Context, param: click.Parameter, value: t.Any) -> t.Any:
        if not value:
            return value
        if (is_flag and bool(value) is True) or (value is not None):
            message = (
                f"'{optstring}' is deprecated and will be removed in a future release."
            )
            if append_message is not None:
                message += f" {append_message}"
            warnings.warn(message, stacklevel=2)

        return value

    return callback
