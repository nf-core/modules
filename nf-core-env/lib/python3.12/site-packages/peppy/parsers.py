import os
from typing import Any, Dict, List

import pandas as pd

from .exceptions import InvalidSampleTableFileException


class TableParser:
    """
    Generic class for sample table parsers.

    Each parser must implement the following methods:
        - parse
        - ...
    """

    def __init__(self, path: str, exts: List[str]) -> None:
        self._path = path
        self._exts = exts
        self._table: pd.DataFrame = None
        self._pandas_kwargs: Dict[str, Any] = {
            "dtype": str,
            "index_col": False,
            "keep_default_na": False,
            "na_values": [""],
        }

    @property
    def path(self) -> str:
        """
        Return the path to the sample table
        """
        return self._path

    @property
    def extensions(self) -> List[str]:
        """
        Return the list of extensions supported by the parser
        """
        return self._exts

    @property
    def table(self) -> pd.DataFrame:
        """
        The parsed table

        returns pandas.DataFrame: the parsed sample table
        """
        return self._table if self._table is not None else self.parse()

    def validate_path(self) -> None:
        """
        Validate the path to the sample table

        Validation includes:
            - check whether extension is supported
            - check whether the file exists
        """
        if not any(self.path.endswith(ext) for ext in self.extensions):
            raise InvalidSampleTableFileException(
                f"Sample table file format not supported: {self.path}"
            )

    def parse(self) -> pd.DataFrame:
        """
        Parse the sample table
        """
        raise NotImplementedError

    def __repr__(self) -> str:
        """
        Return a string representation of the parser
        """
        return f"<{self.__class__.__name__}(path={self.path})>"


class CSVTableParser(TableParser):
    """
    Parser for CSV sample tables
    """

    def __init__(self, path: str) -> None:
        super().__init__(path, ["csv"])

    def parse(self) -> pd.DataFrame:
        """
        Parse the sample table
        """
        self.validate_path()
        self._table = pd.read_csv(self.path, **self._pandas_kwargs)
        self._table = self._table.where(pd.notnull(self._table), None)
        return self.table


class TSVTableParser(TableParser):
    """
    Parser for TSV sample tables
    """

    def __init__(self, path: str) -> None:
        super().__init__(path, ["tsv"])

    def parse(self) -> pd.DataFrame:
        """
        Parse the sample table
        """
        self.validate_path()
        self._table = pd.read_csv(self.path, sep="\t", **self._pandas_kwargs)
        return self.table


class XLSXTableParser(TableParser):
    """
    Parser for MS Excel sample tables
    """

    def __init__(self, path: str) -> None:
        super().__init__(path, ["xlsx"])

    def parse(self) -> pd.DataFrame:
        """
        Parse the sample table
        """
        self.validate_path()
        self._table = pd.read_excel(self.path, **self._pandas_kwargs)
        return self.table


def select_parser(path: str) -> TableParser:
    """
    Select a parser based on the file extension

    :param str path: file path
    :return SampleTableParser: the selected parser
    :raises InvalidSampleTableFileException: if no parser is found for the extension
    """
    parsers_by_ext = parser_by_ext()
    ext = os.path.splitext(path)[1].split(".")[-1]
    if ext in parsers_by_ext:
        return parsers_by_ext[ext]
    raise InvalidSampleTableFileException(
        f"No parser found for extension: {ext}. "
        f"Supported sample table extensions: {list(parsers_by_ext.keys())}",
    )


def parser_by_ext() -> Dict[str, TableParser]:
    """
    Return a dict of parsers indexed by extension

    :return Dict[str, SampleTableParser]: dict of parsers indexed by extension
    """
    parsers_by_ext = {}
    for parser in [cls for cls in TableParser.__subclasses__()]:
        for ext in parser("").extensions:
            parsers_by_ext[ext] = parser
    return parsers_by_ext
