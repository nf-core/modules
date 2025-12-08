from ._version import __version__
from .alias import *

# Origina version
from .yacman import *

# For transition, mostly backwards-compatible version
from .yacman1 import YAMLConfigManager, select_config, load_yaml

# Future version (not backwards-compatible)
from .yacman_future import FutureYAMLConfigManager
from ubiquerg import read_lock, write_lock
