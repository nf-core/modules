import logging
from pathlib import Path

log = logging.getLogger(__name__)


def local_component_structure(self):
    """
    Check that the local modules and subworkflows directories in a pipeline have the correct format:

    .. code-block:: bash

        modules/local/TOOL/SUBTOOL

    Prior to nf-core/tools release 3.1.0 the directory structure allowed top-level `*.nf` files:

    .. code-block:: bash

        modules/local/modules/TOOL_SUBTOOL.nf
    """
    warned_mods = []
    for nf_file in Path(self.wf_path, "modules", "local").glob("*.nf"):
        warned_mods.append(f"{nf_file.name} in modules/local should be moved to a TOOL/SUBTOOL/main.nf structure")
    # If there are modules installed in the wrong location
    passed = []
    if len(warned_mods) == 0:
        passed = ["local modules directory structure is correct 'modules/local/TOOL/SUBTOOL'"]

    warned_swfs = []
    for nf_file in Path(self.wf_path, "subworkflows", "local").glob("*.nf"):
        warned_swfs.append(
            f"{nf_file.name} in subworkflows/local should be moved to a SUBWORKFLOW_NAME/main.nf structure"
        )

    if len(warned_swfs) == 0:
        passed = ["local subworkflows directory structure is correct 'subworkflows/local/TOOL/SUBTOOL'"]

    return {"passed": passed, "warned": warned_mods + warned_swfs, "failed": [], "ignored": []}
