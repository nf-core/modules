"""
Verify that a module has a correct entry in the modules.json file
"""

import logging
from pathlib import Path

import nf_core
import nf_core.modules.lint
import nf_core.modules.modules_repo
import nf_core.modules.modules_utils
from nf_core.modules.modules_utils import NFCoreComponent

log = logging.getLogger(__name__)


def module_version(module_lint_object: "nf_core.modules.lint.ModuleLint", module: "NFCoreComponent"):
    """
    Verifies that the module has a version specified in the ``modules.json`` file

    It checks whether the module has an entry in the ``modules.json`` file
    containing a commit SHA. If that is true, it verifies that there are no
    newer version of the module available.
    """
    assert module_lint_object.modules_json is not None  # mypy
    assert module.repo_url is not None  # mypy
    modules_json_path = Path(module_lint_object.directory, "modules.json")
    # Verify that a git_sha exists in the `modules.json` file for this module
    version = module_lint_object.modules_json.get_module_version(module.component_name, module.repo_url, module.org)
    if version is None:
        module.failed.append(("module_version", "git_sha", "No git_sha entry in `modules.json`", modules_json_path))
        return

    module.git_sha = version
    module.passed.append(("module_version", "git_sha", "Found git_sha entry in `modules.json`", modules_json_path))

    # Check whether a new version is available
    try:
        module.branch = module_lint_object.modules_json.get_component_branch(
            "modules", module.component_name, module.repo_url, module.org
        )
        modules_repo = nf_core.modules.modules_repo.ModulesRepo(remote_url=module.repo_url, branch=module.branch)

        module_git_log = list(modules_repo.get_component_git_log(module.component_name, "modules"))
        if version == module_git_log[0]["git_sha"]:
            module.passed.append(
                ("module_version", "module_version", "Module is the latest version", module.component_dir)
            )
        else:
            module.warned.append(("module_version", "module_version", "New version available", module.component_dir))
    except UserWarning:
        module.warned.append(("module_version", "module_version", "Failed to fetch git log", module.component_dir))
