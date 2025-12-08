from __future__ import annotations

import linecache
import textwrap
import traceback
import typing as t

import click


def format_error_message(err: BaseException) -> str:
    return f"{type(err).__name__}: {err}"


def format_minimal_error(err: BaseException, *, indent: int = 0) -> str:
    lines = [textwrap.indent(str(err), indent * " ")]
    if err.__cause__ is not None:
        lines.append(
            textwrap.indent(format_error_message(err.__cause__), (indent + 2) * " ")
        )

    return "\n".join(lines)


def format_shortened_error(err: BaseException, *, indent: int = 0) -> str:
    lines = []
    lines.append(textwrap.indent(format_error_message(err), indent * " "))
    if err.__traceback__ is not None:
        lineno = err.__traceback__.tb_lineno
        tb_frame = err.__traceback__.tb_frame
        filename = tb_frame.f_code.co_filename
        line = linecache.getline(filename, lineno)
        lines.append((indent + 2) * " " + f'in "{filename}", line {lineno}')
        lines.append((indent + 2) * " " + ">>> " + line.strip())
    return "\n".join(lines)


def format_shortened_trace(caught_err: BaseException) -> str:
    err_stack: list[BaseException] = [caught_err]
    while err_stack[-1].__context__ is not None:
        err_stack.append(err_stack[-1].__context__)  # type: ignore[arg-type]

    parts = [format_shortened_error(caught_err)]
    indent = 0
    for err in err_stack[1:]:
        indent += 2
        parts.append("\n" + indent * " " + "caused by\n")
        parts.append(format_shortened_error(err, indent=indent))
    return "\n".join(parts)


def format_error(
    err: Exception, mode: t.Literal["minimal", "short", "full"] = "short"
) -> str:
    if mode == "minimal":
        return format_minimal_error(err)
    elif mode == "short":
        return format_shortened_trace(err)
    else:
        return "".join(traceback.format_exception(type(err), err, err.__traceback__))


def print_error(
    err: Exception, mode: t.Literal["minimal", "short", "full"] = "short"
) -> None:
    click.echo(format_error(err, mode=mode), err=True)
