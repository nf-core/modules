from pathlib import Path

from ...components.components_differ import ComponentsDiffer
from ...components.nfcore_component import NFCoreComponent


def module_patch(module_lint_obj, module: NFCoreComponent):
    """
    Lint a patch file found in a module

    Checks that the file name is well formed, and that
    the patch can be applied in reverse with the correct result.
    """
    # Check if the module is patched
    if not module.is_patched:
        # Nothing to do
        return

    if not check_patch_valid(module, module.patch_path):
        # Test failed, just exit
        return

    patch_reversible(module_lint_obj, module, module.patch_path)


def check_patch_valid(module, patch_path):
    """
    Checks whether a patch is valid. Looks for lines like
    --- <path>
    +++ <path>
    @@ n,n n,n @@
    and make sure that the come in the right order and that
    the reported paths exists. If the patch file performs
    file creation or deletion we issue a lint warning.

    Args:
        module (NFCoreComponent): The module currently being linted
        patch_path (Path): The absolute path to the patch file.

    Returns:
        (bool): False if any test failed, True otherwise
    """
    with open(patch_path) as fh:
        patch_lines = fh.readlines()

    # Check that the file contains a patch for at least one file
    # and that the file is in the correct directory
    paths_in_patch = []
    passed = True
    it = iter(patch_lines)
    try:
        while True:
            line = next(it)
            if line.startswith("---"):
                frompath = Path(line.split(" ")[1].strip("\n"))
                line = next(it)
                if not line.startswith("+++"):
                    module.failed.append(
                        (
                            "module_patch",
                            "patch_valid",
                            "Patch file invalid. Line starting with '---' should always be followed by line starting with '+++'",
                            patch_path,
                        )
                    )
                    passed = False
                    continue
                topath = Path(line.split(" ")[1].strip("\n"))
                if frompath == Path("/dev/null"):
                    paths_in_patch.append((frompath, ComponentsDiffer.DiffEnum.CREATED))
                elif topath == Path("/dev/null"):
                    paths_in_patch.append((frompath, ComponentsDiffer.DiffEnum.REMOVED))
                elif frompath == topath:
                    paths_in_patch.append((frompath, ComponentsDiffer.DiffEnum.CHANGED))
                else:
                    module.failed.append(
                        (
                            "module_patch",
                            "patch_valid",
                            f"Patch file invaldi. From file '{frompath}' mismatched with to path '{topath}'",
                            patch_path,
                        )
                    )
                    passed = False
                # Check that the next line is hunk
                line = next(it)
                if not line.startswith("@@"):
                    module.failed.append(
                        (
                            "module_patch",
                            "patch_valid",
                            "Patch file invalid. File declarations should be followed by hunk",
                            patch_path,
                        )
                    )
                    passed = False
    except StopIteration:
        pass

    if not passed:
        return False

    if len(paths_in_patch) == 0:
        module.failed.append(("module_patch", "patch_valid", "Patch file invalid. Found no patches", patch_path))
        return False

    # Go through the files and check that they exist
    # Warn about any created or removed files
    passed = True
    for path, diff_status in paths_in_patch:
        if diff_status == ComponentsDiffer.DiffEnum.CHANGED:
            if not Path(module.base_dir, path).exists():
                module.failed.append(
                    (
                        "module_patch",
                        "patch_valid",
                        f"Patch file invalid. Path '{path}' does not exist but is reported in patch file.",
                        patch_path,
                    )
                )
                passed = False
                continue
        elif diff_status == ComponentsDiffer.DiffEnum.CREATED:
            if not Path(module.base_dir, path).exists():
                module.failed.append(
                    (
                        "module_patch",
                        "patch_valid",
                        f"Patch file invalid. Path '{path}' does not exist but is reported in patch file.",
                        patch_path,
                    )
                )
                passed = False
                continue
            module.warned.append(
                ("module_patch", "patch", f"Patch file performs file creation of {path}. This is discouraged."),
                patch_path,
            )
        elif diff_status == ComponentsDiffer.DiffEnum.REMOVED:
            if Path(module.base_dir, path).exists():
                module.failed.append(
                    (
                        "module_patch",
                        "patch_valid",
                        f"Patch file invalid. Path '{path}' is reported as deleted but exists.",
                        patch_path,
                    )
                )
                passed = False
                continue
            module.warned.append(
                (
                    "module_patch",
                    "patch",
                    f"Patch file performs file deletion of {path}. This is discouraged.",
                    patch_path,
                )
            )
        if passed:
            module.passed.append(("module_patch", "patch_valid", "Patch file is valid", patch_path))
        return passed


def patch_reversible(module_lint_object, module, patch_path):
    """
    Try applying a patch in reverse to see if it is up to date

    Args:
        module (NFCoreComponent): The module currently being linted
        patch_path (Path): The absolute path to the patch file.

    Returns:
        (bool): False if any test failed, True otherwise
    """
    try:
        ComponentsDiffer.try_apply_patch(
            module.component_type,
            module.component_name,
            module_lint_object.modules_repo.repo_path,
            patch_path,
            Path(module.component_dir).relative_to(module.base_dir),
            reverse=True,
        )
    except LookupError:
        # Patch failed. Save the patch file by moving to the install dir
        module.failed.append(("module_patch", "patch_reversible", "Patch file is outdated or edited", patch_path))
        return False

    module.passed.append(("module_patch", "patch_reversible", "Patch agrees with module files", patch_path))
    return True
