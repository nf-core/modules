from __future__ import annotations

import io
import typing as t

from check_jsonschema.cli.param_types import CustomLazyFile

from .parsers import ParseError, ParserSet
from .transforms import Transform


class InstanceLoader:
    def __init__(
        self,
        files: t.Sequence[t.IO[bytes] | CustomLazyFile],
        default_filetype: str = "json",
        force_filetype: str | None = None,
        data_transform: Transform | None = None,
    ) -> None:
        self._files = files
        self._default_filetype = default_filetype
        self._force_filetype = force_filetype
        self._data_transform = (
            data_transform if data_transform is not None else Transform()
        )

        self._parsers = ParserSet(
            modify_yaml_implementation=self._data_transform.modify_yaml_implementation
        )

    def iter_files(self) -> t.Iterator[tuple[str, ParseError | t.Any]]:
        for file in self._files:
            if hasattr(file, "name"):
                name = file.name
            # allowing for BytesIO to be special-cased here is useful for
            # simpler test setup, since this is what tests will pass and we naturally
            # support it here
            elif isinstance(file, io.BytesIO) or file.fileno() == 0:
                name = "<stdin>"
            else:
                raise ValueError(f"File {file} has no name attribute")

            try:
                if isinstance(file, CustomLazyFile):
                    stream: t.IO[bytes] = t.cast(t.IO[bytes], file.open())
                else:
                    stream = file

                try:
                    data: t.Any = self._parsers.parse_data_with_path(
                        stream, name, self._default_filetype, self._force_filetype
                    )
                except ParseError as err:
                    data = err
                else:
                    data = self._data_transform(data)
            finally:
                file.close()
            yield (name, data)
