import difflib
import enum
import json
import logging
import os
from pathlib import Path

from rich import box
from rich.console import Console, Group, RenderableType
from rich.panel import Panel
from rich.syntax import Syntax

import nf_core.utils

log = logging.getLogger(__name__)


class ComponentsDiffer:
    """
    Static class that provides functionality for computing diffs between
    different instances of a module or subworkflow
    """

    class DiffEnum(enum.Enum):
        """
        Enumeration for keeping track of
        the diff status of a pair of files
        """

        UNCHANGED = enum.auto()
        CHANGED = enum.auto()
        CREATED = enum.auto()
        REMOVED = enum.auto()

    @staticmethod
    def get_component_diffs(from_dir, to_dir, for_git=True, dsp_from_dir=None, dsp_to_dir=None):
        """
        Compute the diff between the current component version
        and the new version.

        Args:
            from_dir (strOrPath): The folder containing the old component files
            to_dir (strOrPath): The folder containing the new component files
            path_in_diff (strOrPath): The directory displayed containing the component
                                      file in the diff. Added so that temporary dirs
                                      are not shown
            for_git (bool): indicates whether the diff file is to be
                            compatible with `git apply`. If true it
                            adds a/ and b/ prefixes to the file paths
            dsp_from_dir (str | Path): The from directory to display in the diff
            dsp_to_dir (str | Path): The to directory to display in the diff

        Returns:
            dict[str, (ComponentsDiffer.DiffEnum, str)]: A dictionary containing
            the diff type and the diff string (empty if no diff)
        """
        if for_git:
            dsp_from_dir = Path("a", dsp_from_dir)
            dsp_to_dir = Path("b", dsp_to_dir)

        diffs = {}
        # Get all unique filenames in the two folders.
        # `dict.fromkeys()` is used instead of `set()` to preserve order
        files = dict.fromkeys(
            [Path(dirpath, file).relative_to(to_dir) for dirpath, _, files in os.walk(to_dir) for file in files]
        )
        files.update(
            dict.fromkeys(
                [Path(dirpath, file).relative_to(from_dir) for dirpath, _, files in os.walk(from_dir) for file in files]
            )
        )
        files = list(files)

        # Loop through all the component files and compute their diffs if needed
        for file in files:
            temp_path = Path(to_dir, file)
            curr_path = Path(from_dir, file)
            if temp_path.exists() and curr_path.exists() and temp_path.is_file():
                with open(temp_path) as fh:
                    new_lines = fh.readlines()
                with open(curr_path) as fh:
                    old_lines = fh.readlines()

                if new_lines == old_lines:
                    # The files are identical
                    diffs[file] = (ComponentsDiffer.DiffEnum.UNCHANGED, ())
                else:
                    # Compute the diff
                    diff = difflib.unified_diff(
                        old_lines,
                        new_lines,
                        fromfile=str(Path(dsp_from_dir, file)),
                        tofile=str(Path(dsp_to_dir, file)),
                    )
                    diffs[file] = (ComponentsDiffer.DiffEnum.CHANGED, diff)

            elif temp_path.exists():
                with open(temp_path) as fh:
                    new_lines = fh.readlines()
                # The file was created
                # Show file against /dev/null
                diff = difflib.unified_diff(
                    [],
                    new_lines,
                    fromfile=str(Path("/dev", "null")),
                    tofile=str(Path(dsp_to_dir, file)),
                )
                diffs[file] = (ComponentsDiffer.DiffEnum.CREATED, diff)

            elif curr_path.exists():
                # The file was removed
                # Show file against /dev/null
                with open(curr_path) as fh:
                    old_lines = fh.readlines()
                diff = difflib.unified_diff(
                    old_lines,
                    [],
                    fromfile=str(Path(dsp_from_dir, file)),
                    tofile=str(Path("/dev", "null")),
                )
                diffs[file] = (ComponentsDiffer.DiffEnum.REMOVED, diff)

        return diffs

    @staticmethod
    def write_diff_file(
        diff_path,
        component,
        repo_path,
        from_dir,
        to_dir,
        current_version=None,
        new_version=None,
        file_action="a",
        for_git=True,
        dsp_from_dir=None,
        dsp_to_dir=None,
        limit_output=False,
    ):
        """
        Writes the diffs of a component to the diff file.

        Args:
            diff_path (str | Path): The path to the file that should be appended
            component (str): The component name
            repo_path (str): The name of the repo where the component resides
            from_dir (str | Path): The directory containing the old component files
            to_dir (str | Path): The directory containing the new component files
            diffs (dict[str, (ComponentsDiffer.DiffEnum, str)]): A dictionary containing
                                                              the type of change and
                                                              the diff (if any)
            current_version (str): The installed version of the component
            new_version (str): The version of the component the diff is computed against
            for_git (bool): indicates whether the diff file is to be
                            compatible with `git apply`. If true it
                            adds a/ and b/ prefixes to the file paths
            dsp_from_dir (str | Path): The 'from' directory displayed in the diff
            dsp_to_dir (str | Path): The 'to' directory displayed in the diff
            limit_output (bool): If true, don't write the diff for files other than main.nf
        """
        if dsp_from_dir is None:
            dsp_from_dir = from_dir
        if dsp_to_dir is None:
            dsp_to_dir = to_dir

        diffs = ComponentsDiffer.get_component_diffs(from_dir, to_dir, for_git, dsp_from_dir, dsp_to_dir)
        if all(diff_status == ComponentsDiffer.DiffEnum.UNCHANGED for _, (diff_status, _) in diffs.items()):
            raise UserWarning("Component is unchanged")
        log.debug(f"Writing diff of '{component}' to '{diff_path}'")
        with open(diff_path, file_action) as fh:
            if current_version is not None and new_version is not None:
                fh.write(
                    f"Changes in component '{Path(repo_path, component)}' between"
                    f" ({current_version}) and"
                    f" ({new_version})\n"
                )
            else:
                fh.write(f"Changes in component '{Path(repo_path, component)}'\n")

            for file, (diff_status, diff) in diffs.items():
                if diff_status == ComponentsDiffer.DiffEnum.UNCHANGED:
                    # The files are identical
                    fh.write(f"'{Path(dsp_from_dir, file)}' is unchanged\n")
                elif diff_status == ComponentsDiffer.DiffEnum.CREATED:
                    # The file was created between the commits
                    fh.write(f"'{Path(dsp_from_dir, file)}' was created\n")
                elif diff_status == ComponentsDiffer.DiffEnum.REMOVED:
                    # The file was removed between the commits
                    fh.write(f"'{Path(dsp_from_dir, file)}' was removed\n")
                elif limit_output and not file.suffix == ".nf":
                    # Skip printing the diff for files other than main.nf
                    fh.write(f"Changes in '{Path(component, file)}' but not shown\n")
                else:
                    # The file has changed write the diff lines to the file
                    fh.write(f"Changes in '{Path(component, file)}':\n")
                    for line in diff:
                        fh.write(line)
                    fh.write("\n")

            fh.write("*" * 60 + "\n")

    @staticmethod
    def append_modules_json_diff(diff_path, old_modules_json, new_modules_json, modules_json_path, for_git=True):
        """
        Compare the new modules.json and builds a diff

        Args:
            diff_fn (str): The diff file to be appended
            old_modules_json (nested dict): The old modules.json
            new_modules_json (nested dict): The new modules.json
            modules_json_path (str): The path to the modules.json
            for_git (bool): indicates whether the diff file is to be
                            compatible with `git apply`. If true it
                            adds a/ and b/ prefixes to the file paths
        """
        fromfile = modules_json_path
        tofile = modules_json_path
        if for_git:
            fromfile = Path("a", fromfile)
            tofile = Path("b", tofile)

        modules_json_diff = difflib.unified_diff(
            json.dumps(old_modules_json, indent=4).splitlines(keepends=True),
            json.dumps(new_modules_json, indent=4).splitlines(keepends=True),
            fromfile=str(fromfile),
            tofile=str(tofile),
        )

        # Save diff for modules.json to file
        with open(diff_path, "a") as fh:
            fh.write("Changes in './modules.json'\n")
            for line in modules_json_diff:
                fh.write(line)
            fh.write("*" * 60 + "\n")

    @staticmethod
    def print_diff(
        component,
        repo_path,
        from_dir,
        to_dir,
        current_version=None,
        new_version=None,
        dsp_from_dir=None,
        dsp_to_dir=None,
        limit_output=False,
    ):
        """
        Prints the diffs between two component versions to the terminal

        Args:
            component (str): The component name
            repo_path (str): The name of the repo where the component resides
            from_dir (str | Path): The directory containing the old component files
            to_dir (str | Path): The directory containing the new component files
            current_version (str): The installed version of the component
            new_version (str): The version of the component the diff is computed against
            dsp_from_dir (str | Path): The 'from' directory displayed in the diff
            dsp_to_dir (str | Path): The 'to' directory displayed in the diff
            limit_output (bool): If true, don't print the diff for files other than main.nf
        """
        if dsp_from_dir is None:
            dsp_from_dir = from_dir
        if dsp_to_dir is None:
            dsp_to_dir = to_dir

        diffs = ComponentsDiffer.get_component_diffs(
            from_dir, to_dir, for_git=False, dsp_from_dir=dsp_from_dir, dsp_to_dir=dsp_to_dir
        )
        console = Console(force_terminal=nf_core.utils.rich_force_colors())
        if current_version is not None and new_version is not None:
            log.info(
                f"Changes in component '{Path(repo_path, component)}' between ({current_version}) and ({new_version})"
            )
        else:
            log.info(f"Changes in component '{Path(repo_path, component)}'")

        panel_group: list[RenderableType] = []
        for file, (diff_status, diff) in diffs.items():
            if diff_status == ComponentsDiffer.DiffEnum.UNCHANGED:
                # The files are identical
                log.info(f"'{Path(dsp_from_dir, file)}' is unchanged")
            elif diff_status == ComponentsDiffer.DiffEnum.CREATED:
                # The file was created between the commits
                log.info(f"'{Path(dsp_from_dir, file)}' was created")
            elif diff_status == ComponentsDiffer.DiffEnum.REMOVED:
                # The file was removed between the commits
                log.info(f"'{Path(dsp_from_dir, file)}' was removed")
            elif limit_output and not file.suffix == ".nf":
                # Skip printing the diff for files other than main.nf
                log.info(f"Changes in '{Path(component, file)}' but not shown")
            else:
                # The file has changed
                log.info(f"Changes in '{Path(component, file)}':")
                # Pretty print the diff using the pygments diff lexer
                syntax = Syntax("".join(diff), "diff", theme="ansi_dark", line_numbers=True)
                panel_group.append(Panel(syntax, title=str(file), title_align="left", padding=0))
        console.print(
            Panel(
                Group(*panel_group),
                title=f"[white]{str(component)}[/white]",
                title_align="left",
                padding=0,
                border_style="blue",
                box=box.HEAVY,
            )
        )

    @staticmethod
    def per_file_patch(patch_fn: str | Path) -> dict[str, list[str]]:
        """
        Splits a patch file for several files into one patch per file.

        Args:
            patch_fn (str | Path): The path to the patch file

        Returns:
            dict[str, str]: A dictionary indexed by the filenames with the
                            file patches as values
        """
        with open(patch_fn) as fh:
            lines = fh.readlines()

        patches = {}
        i = 0
        patch_lines: list[str] = []
        key = "preamble"
        while i < len(lines):
            line = lines[i]
            if line.startswith("---"):
                # New file found: add the old lines to the dictionary
                # and determine the new filename
                patches[key] = patch_lines
                _, frompath = line.split(" ")
                frompath = frompath.strip()
                patch_lines = [line]
                i += 1
                line = lines[i]
                if not line.startswith("+++ "):
                    raise LookupError("Missing to-file in patch file")
                _, topath = line.split(" ")
                topath = topath.strip()
                patch_lines.append(line)

                if frompath == topath:
                    key = frompath
                elif frompath == "/dev/null":
                    key = topath
                else:
                    key = frompath
            else:
                patch_lines.append(line)
            i += 1
        patches[key] = patch_lines

        # Remove the 'preamble' key (entry contains no useful information)
        patches.pop("preamble")
        return patches

    @staticmethod
    def get_new_and_old_lines(patch):
        """
        Parse a patch for a file, and return the contents
        of the modified parts for both the old and new versions

        Args:
            patch (str): The patch in unified diff format

        Returns:
            ([[str]], [[str]]): Lists of old and new lines for each hunk
                                (modified part the file)
        """
        old_lines = []
        new_lines = []
        old_partial = []
        new_partial = []
        # First two lines indicate the file names, the third line is the first hunk
        for line in patch[3:]:
            if line.startswith("@"):
                old_lines.append(old_partial)
                new_lines.append(new_partial)
                old_partial = []
                new_partial = []
            elif line.startswith(" "):
                # Line belongs to both files
                line = line[1:]
                old_partial.append(line)
                new_partial.append(line)
            elif line.startswith("+"):
                # Line only belongs to the new file
                line = line[1:]
                new_partial.append(line)
            elif line.startswith("-"):
                # Line only belongs to the old file
                line = line[1:]
                old_partial.append(line)
        old_lines.append(old_partial)
        new_lines.append(new_partial)
        return old_lines, new_lines

    @staticmethod
    def try_apply_single_patch(file_lines, patch, reverse=False):
        """
        Tries to apply a patch to a modified file. Since the line numbers in
        the patch do not agree if the file is modified, the old and new
        lines in the patch are reconstructed and then we look for the old lines
        in the modified file. If all hunk in the patch are found in the new file
        it is updated with the new lines from the patch file.

        Args:
            file_lines ([str]): The lines of the file to be patched
            patch (str | Path): (Outdated) patch for the file
            reverse (bool): Apply the patch in reverse

        Returns:
            [str]: The patched lines of the file

        Raises:
            LookupError: If it fails to find the old lines from the patch in
                         the file.
        """
        org_lines, patch_lines = ComponentsDiffer.get_new_and_old_lines(patch)
        if reverse:
            patch_lines, org_lines = org_lines, patch_lines

        # The patches are sorted by their order of occurrence in the original
        # file. Loop through the new file and try to find the new indices of
        # these lines. We know they are non overlapping, and thus only need to
        # look at the file once
        p = len(org_lines)
        patch_indices = [None] * p
        i = 0
        j = 0
        n = len(file_lines)
        while i < n and j < p:
            m = len(org_lines[j])
            while i < n:
                if org_lines[j] == file_lines[i : i + m]:
                    patch_indices[j] = (i, i + m)
                    j += 1
                    break
                else:
                    i += 1

        if j != len(org_lines):
            # We did not find all diffs before we ran out of file.
            raise LookupError("Failed to find lines where patch should be applied")

        # Apply the patch to new lines by substituting
        # the original lines with the patch lines
        patched_new_lines = file_lines[: patch_indices[0][0]]
        for i in range(len(patch_indices) - 1):
            # Add the patch lines
            patched_new_lines.extend(patch_lines[i])
            # Fill the spaces between the patches
            patched_new_lines.extend(file_lines[patch_indices[i][1] : patch_indices[i + 1][0]])

        # Add the remaining part of the new file
        patched_new_lines.extend(patch_lines[-1])
        patched_new_lines.extend(file_lines[patch_indices[-1][1] :])

        return patched_new_lines

    @staticmethod
    def try_apply_patch(
        component_type: str,
        component: str,
        repo_path: str | Path,
        patch_path: str | Path,
        component_dir: Path,
        reverse: bool = False,
    ) -> dict[str, list[str]]:
        """
        Try applying a full patch file to a module or subworkflow

        Args:
            component_type (str): The type of component (modules or subworkflows)
            component (str): Name of the module or subworkflow
            repo_path (str): Name of the repository where the component resides
            patch_path (str): The absolute path to the patch file to be applied
            component_dir (Path): The directory containing the component
            reverse (bool): Apply the patch in reverse

        Returns:
            dict[str, str]: A dictionary with file paths (relative to the pipeline dir)
                            as keys and the patched file contents as values

        Raises:
            LookupError: If the patch application fails in a file
        """
        component_relpath = Path(component_type, repo_path, component)
        patches = ComponentsDiffer.per_file_patch(patch_path)
        new_files = {}
        for file, patch in patches.items():
            log.debug(f"Applying patch to {file}")
            fn = Path(file).relative_to(component_relpath)
            file_path = component_dir / fn
            try:
                with open(file_path) as fh:
                    file_lines = fh.readlines()
            except FileNotFoundError:
                # The file was added with the patch
                file_lines = [""]
            patched_new_lines = ComponentsDiffer.try_apply_single_patch(file_lines, patch, reverse=reverse)
            new_files[str(fn)] = patched_new_lines
        return new_files
