import logging
from pathlib import Path

from nf_core.components.list import ComponentList

log = logging.getLogger(__name__)


class ModuleList(ComponentList):
    def __init__(
        self,
        pipeline_dir: str | Path = ".",
        remote: bool = True,
        remote_url: str | None = None,
        branch: str | None = None,
        no_pull: bool = False,
    ):
        super().__init__("modules", pipeline_dir, remote, remote_url, branch, no_pull)
