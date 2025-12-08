import logging
import os
from pathlib import Path

log = logging.getLogger(__name__)


def modules_structure(self):
    """
    Check that the structure of the modules directory in a pipeline is the correct one:

    .. code-block:: bash

        modules/nf-core/TOOL/SUBTOOL

    Prior to nf-core/tools release 2.6 the directory structure had an additional level of nesting:

    .. code-block:: bash

        modules/nf-core/modules/TOOL/SUBTOOL
    """
    wrong_location_modules = []
    for directory, _, files in os.walk(Path(self.wf_path, "modules")):
        if "main.nf" in files:
            module_path = Path(directory).relative_to(Path(self.wf_path, "modules"))
            parts = module_path.parts
            # Check that there are modules installed directly under the 'modules' directory
            if parts[1] == "modules":
                wrong_location_modules.append(module_path)
    # If there are modules installed in the wrong location
    failed = []
    passed = []
    if len(wrong_location_modules) > 0:
        failed = ["modules directory structure is outdated. Should be 'modules/nf-core/TOOL/SUBTOOL'"]
    else:
        passed = ["modules directory structure is correct 'modules/nf-core/TOOL/SUBTOOL'"]
    return {"passed": passed, "warned": [], "failed": failed, "ignored": []}
