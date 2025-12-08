import logging

from nf_core.components.create import ComponentCreate

log = logging.getLogger(__name__)


class ModuleCreate(ComponentCreate):
    def __init__(
        self,
        pipeline_dir,
        component="",
        author=None,
        process_label=None,
        has_meta=None,
        force=False,
        conda_name=None,
        conda_version=None,
        empty_template=False,
        migrate_pytest=False,
    ):
        super().__init__(
            "modules",
            pipeline_dir,
            component,
            author,
            process_label,
            has_meta,
            force,
            conda_name,
            conda_version,
            empty_template,
            migrate_pytest,
        )
