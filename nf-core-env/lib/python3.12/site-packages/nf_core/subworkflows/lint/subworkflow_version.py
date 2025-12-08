"""
Verify that a subworkflow has a correct entry in the modules.json file
"""

import logging
from pathlib import Path

import nf_core
import nf_core.modules.modules_repo
import nf_core.modules.modules_utils

log = logging.getLogger(__name__)


def subworkflow_version(subworkflow_lint_object, subworkflow):
    """
    Verifies that the subworkflow has a version specified in the ``modules.json`` file

    It checks whether the subworkflow has an entry in the ``modules.json`` file
    containing a commit SHA. If that is true, it verifies that there are no
    newer version of the subworkflow available.
    """

    modules_json_path = Path(subworkflow_lint_object.directory, "modules.json")
    # Verify that a git_sha exists in the `modules.json` file for this module
    version = subworkflow_lint_object.modules_json.get_subworkflow_version(
        subworkflow.component_name, subworkflow.repo_url, subworkflow.org
    )
    if version is None:
        subworkflow.failed.append(
            ("subworkflow_version", "git_sha", "No git_sha entry in `modules.json`", modules_json_path)
        )
        return

    subworkflow.git_sha = version
    subworkflow.passed.append(
        ("subworkflow_version", "git_sha", "Found git_sha entry in `modules.json`", modules_json_path)
    )

    # Check whether a new version is available
    try:
        subworkflow.branch = subworkflow_lint_object.modules_json.get_component_branch(
            "subworkflows", subworkflow.component_name, subworkflow.repo_url, subworkflow.org
        )
        modules_repo = nf_core.modules.modules_repo.ModulesRepo(
            remote_url=subworkflow.repo_url, branch=subworkflow.branch
        )

        subworkflow_git_log = modules_repo.get_component_git_log(subworkflow.component_name, "subworkflows")
        if version == next(subworkflow_git_log)["git_sha"]:
            subworkflow.passed.append(
                (
                    "subworkflow_version",
                    "subworkflow_version",
                    "Subworkflow is in the latest version",
                    subworkflow.component_dir,
                )
            )
        else:
            subworkflow.warned.append(
                ("subworkflow_version", "subworkflow_version", "New version available", subworkflow.component_dir)
            )
    except UserWarning:
        subworkflow.warned.append(
            ("subworkflow_version", "subworkflow_version", "Failed to fetch git log", subworkflow.component_dir)
        )
