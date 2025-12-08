import logging

from nf_core.components.info import ComponentInfo

log = logging.getLogger(__name__)


class SubworkflowInfo(ComponentInfo):
    def __init__(
        self,
        pipeline_dir,
        component_name,
        remote_url=None,
        branch=None,
        no_pull=False,
    ):
        super().__init__("subworkflows", pipeline_dir, component_name, remote_url, branch, no_pull)
