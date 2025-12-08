""" Package-scope definitions """

from ._att_map_like import AttMapLike
from ._version import __version__
from .attmap import AttMap
from .attmap_echo import *
from .helpers import *
from .ordattmap import OrdAttMap
from .pathex_attmap import PathExAttMap

AttributeDict = AttMap
AttributeDictEcho = AttMapEcho

__all__ = [
    "AttMapLike",
    "AttMap",
    "AttMapEcho",
    "AttributeDict",
    "AttributeDictEcho",
    "EchoAttMap",
    "OrdAttMap",
    "PathExAttMap",
    "get_data_lines",
]
__aliases__ = {
    "AttMap": ["AttributeDict"],
    "AttMapEcho": ["AttributeDictEcho"],
    "EchoAttMap": ["AttMapEcho"],
}
