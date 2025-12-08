"""
Check whether the content of a module has changed compared to the original repository
"""

import shutil
import tempfile
from pathlib import Path

import nf_core.modules.modules_repo
from nf_core.components.components_differ import ComponentsDiffer


def module_changes(module_lint_object, module):
    """
    Checks whether installed nf-core modules have changed compared to the
    original repository

    Downloads the ``main.nf`` and ``meta.yml`` files for every module
    and compares them to the local copies

    If the module has a commit SHA entry in the ``modules.json``, the file content is
    compared against the files in the remote at this SHA.

    Only runs when linting a pipeline, not the modules repository
    """
    if module.is_patched:
        # If the module is patched, we need to apply
        # the patch in reverse before comparing with the remote
        tempdir_parent = Path(tempfile.mkdtemp())
        tempdir = tempdir_parent / "tmp_module_dir"
        shutil.copytree(module.component_dir, tempdir)
        try:
            new_lines = ComponentsDiffer.try_apply_patch(
                module.component_type,
                module.component_name,
                module.org,
                module.patch_path,
                tempdir,
                reverse=True,
            )
            for file, lines in new_lines.items():
                with open(tempdir / file, "w") as fh:
                    fh.writelines(lines)
        except LookupError:
            # This error is already reported by module_patch, so just return
            return
    else:
        tempdir = module.component_dir
    module.branch = module_lint_object.modules_json.get_component_branch(
        "modules", module.component_name, module.repo_url, module.org
    )
    modules_repo = nf_core.modules.modules_repo.ModulesRepo(remote_url=module.repo_url, branch=module.branch)

    for f, same in modules_repo.component_files_identical(
        module.component_name, tempdir, module.git_sha, "modules"
    ).items():
        if same:
            module.passed.append(
                (
                    "module_changes",
                    "check_local_copy",
                    "Local copy of module up to date",
                    f"{Path(module.component_dir, f)}",
                )
            )
        else:
            module.failed.append(
                (
                    "module_changes",
                    "check_local_copy",
                    "Local copy of module does not match remote",
                    f"{Path(module.component_dir, f)}",
                )
            )
