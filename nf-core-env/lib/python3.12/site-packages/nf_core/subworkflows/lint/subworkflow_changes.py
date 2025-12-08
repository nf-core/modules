"""
Check whether the content of a subworkflow has changed compared to the original repository
"""

import shutil
import tempfile
from pathlib import Path

import nf_core.modules.modules_repo
from nf_core.components.components_differ import ComponentsDiffer


def subworkflow_changes(subworkflow_lint_object, subworkflow):
    """
    Checks whether installed nf-core subworkflow have changed compared to the
    original repository

    Downloads the ``main.nf`` and ``meta.yml`` files for every subworkflow
    and compares them to the local copies

    If the subworkflow has a commit SHA entry in the ``modules.json``, the file content is
    compared against the files in the remote at this SHA.

    Only runs when linting a pipeline, not the modules repository
    """
    if subworkflow.is_patched:
        # If the subworkflow is patched, we need to apply
        # the patch in reverse before comparing with the remote
        tempdir_parent = Path(tempfile.mkdtemp())
        tempdir = tempdir_parent / "tmp_subworkflow_dir"
        shutil.copytree(subworkflow.component_dir, tempdir)
        try:
            new_lines = ComponentsDiffer.try_apply_patch(
                subworkflow.component_type,
                subworkflow.component_name,
                subworkflow.org,
                subworkflow.patch_path,
                tempdir,
                reverse=True,
            )
            for file, lines in new_lines.items():
                with open(tempdir / file, "w") as fh:
                    fh.writelines(lines)
            subworkflow.passed.append(
                (
                    "subworkflow_changes",
                    "subworkflow_patch",
                    "Subworkflow patch can be cleanly applied",
                    f"{subworkflow.component_dir}",
                )
            )
        except LookupError:
            subworkflow.failed.append(
                (
                    "subworkflow_changes",
                    "subworkflow_patch",
                    "Subworkflow patch cannot be cleanly applied",
                    f"{subworkflow.component_dir}",
                )
            )
            return
    else:
        tempdir = subworkflow.component_dir
    subworkflow.branch = subworkflow_lint_object.modules_json.get_component_branch(
        "subworkflows", subworkflow.component_name, subworkflow.repo_url, subworkflow.org
    )
    modules_repo = nf_core.modules.modules_repo.ModulesRepo(remote_url=subworkflow.repo_url, branch=subworkflow.branch)

    for f, same in modules_repo.component_files_identical(
        subworkflow.component_name, tempdir, subworkflow.git_sha, "subworkflows"
    ).items():
        if same:
            subworkflow.passed.append(
                (
                    "subworkflow_changes",
                    "check_local_copy",
                    "Local copy of subworkflow up to date",
                    f"{Path(subworkflow.component_dir, f)}",
                )
            )
        else:
            subworkflow.failed.append(
                (
                    "subworkflow_changes",
                    "check_local_copy",
                    "Local copy of subworkflow does not match remote",
                    f"{Path(subworkflow.component_dir, f)}",
                )
            )
