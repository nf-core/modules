# Project configuration, particularly for logging.

import logmuse

from ._version import __version__
from .const import PKG_NAME
from .exceptions import PipestatError
from .pipestat import PipestatBoss, PipestatManager, ProjectPipestatManager, SamplePipestatManager

__all__ = [
    "PipestatError",
    "SamplePipestatManager",
    "ProjectPipestatManager",
    "PipestatBoss",
    "PipestatManager",
    "__version__",
]

logmuse.init_logger(PKG_NAME)
