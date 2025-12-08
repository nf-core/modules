import logging

from nf_core.components.remove import ComponentRemove

log = logging.getLogger(__name__)


class ModuleRemove(ComponentRemove):
    def __init__(self, pipeline_dir, remote_url=None, branch=None, no_pull=False):
        super().__init__("modules", pipeline_dir, remote_url=remote_url, branch=branch, no_pull=no_pull)
