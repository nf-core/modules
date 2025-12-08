import logging

from nf_core.components.patch import ComponentPatch

log = logging.getLogger(__name__)


class ModulePatch(ComponentPatch):
    def __init__(self, pipeline_dir, remote_url=None, branch=None, no_pull=False, installed_by=None):
        super().__init__(pipeline_dir, "modules", remote_url, branch, no_pull, installed_by)
