import logging
import os

log = logging.getLogger(__name__)


def module_deprecations(_, module):
    """
    Check that the modules are up to the latest nf-core standard
    """
    module.wf_path = module.component_dir
    if "functions.nf" in os.listdir(module.component_dir):
        module.failed.append(
            (
                "module_deprecations",
                "module_deprecations",
                "Deprecated file `functions.nf` found. No longer required for the latest nf-core/modules syntax!",
                module.component_dir,
            )
        )
    else:
        module.passed.append(
            (
                "module_deprecations",
                "module_deprecations",
                "No deprecated file `functions.nf` found",
                module.component_dir,
            )
        )
