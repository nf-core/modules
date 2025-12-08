import json
import logging
import subprocess
from pathlib import Path

import rich
import yaml
from rich.console import Console
from rich.table import Table

import nf_core.utils
from nf_core.utils import plural_s as _s

log = logging.getLogger(__name__)

# Create a console used by all lint tests
console = Console(force_terminal=nf_core.utils.rich_force_colors())


def print_joint_summary(lint_obj, module_lint_obj, subworkflow_lint_obj):
    """Print a joint summary of the general pipe lint tests and the module and subworkflow lint tests"""
    swf_passed = 0
    swf_warned = 0
    swf_failed = 0
    module_passed = 0
    module_warned = 0
    module_failed = 0
    if subworkflow_lint_obj is not None:
        swf_passed = len(subworkflow_lint_obj.passed)
        swf_warned = len(subworkflow_lint_obj.warned)
        swf_failed = len(subworkflow_lint_obj.failed)
    if module_lint_obj is not None:
        module_passed = len(module_lint_obj.passed)
        module_warned = len(module_lint_obj.warned)
        module_failed = len(module_lint_obj.failed)
    nbr_passed = len(lint_obj.passed) + module_passed + swf_passed
    nbr_ignored = len(lint_obj.ignored)
    nbr_fixed = len(lint_obj.fixed)
    nbr_warned = len(lint_obj.warned) + module_warned + swf_warned
    nbr_failed = len(lint_obj.failed) + module_failed + swf_failed

    summary_colour = "red" if nbr_failed > 0 else "green"
    table = Table(box=rich.box.ROUNDED, style=summary_colour)
    table.add_column("LINT RESULTS SUMMARY", no_wrap=True)
    table.add_row(rf"[green][✔] {nbr_passed:>3} Test{_s(nbr_passed)} Passed")
    if nbr_fixed:
        table.add_row(rf"[bright blue][?] {nbr_fixed:>3} Test{_s(nbr_fixed)} Fixed")
    table.add_row(rf"[grey58][?] {nbr_ignored:>3} Test{_s(nbr_ignored)} Ignored")
    table.add_row(rf"[yellow][!] {nbr_warned:>3} Test Warning{_s(nbr_warned)}")
    table.add_row(rf"[red][✗] {nbr_failed:>3} Test{_s(nbr_failed)} Failed")
    console.print(table)


def print_fixes(lint_obj):
    """Prints available and applied fixes"""

    if lint_obj.could_fix:
        fix_flags = "".join([f" --fix {fix}" for fix in lint_obj.could_fix])
        wf_dir = "" if lint_obj.wf_path == "." else f"--dir {lint_obj.wf_path}"
        fix_cmd = f"nf-core pipelines lint {wf_dir} {fix_flags}"
        console.print(
            "\nTip: Some of these linting errors can automatically be resolved with the following command:\n\n"
            f"[blue]    {fix_cmd}\n"
        )
    if len(lint_obj.fix):
        console.print(
            "Automatic fixes applied. "
            "Please check with 'git diff' and revert any changes you do not want with 'git checkout <file>'."
        )


def check_git_repo() -> bool:
    """Check if the current directory is a git repository."""
    try:
        subprocess.check_output(["git", "rev-parse", "--is-inside-work-tree"])
        return True
    except subprocess.CalledProcessError:
        return False


def run_prettier_on_file(file: Path | str | list[str]) -> None:
    """Run the pre-commit hook prettier on a file.

    Args:
        file (Path | str): A file identifier as a string or pathlib.Path.

    Warns:
        If Prettier is not installed, a warning is logged.
    """

    is_git = check_git_repo()

    nf_core_pre_commit_config = Path(nf_core.__file__).parent / ".pre-commit-prettier-config.yaml"
    args = ["pre-commit", "run", "--config", str(nf_core_pre_commit_config), "prettier"]
    if isinstance(file, list):
        args.extend(["--files", *file])
    else:
        args.extend(["--files", str(file)])

    if is_git:
        try:
            proc = subprocess.run(args, capture_output=True, check=True)
            log.debug(f"{proc.stdout.decode()}")
        except subprocess.CalledProcessError as e:
            if ": SyntaxError: " in e.stdout.decode():
                log.critical(f"Can't format {file} because it has a syntax error.\n{e.stdout.decode()}")
            elif "files were modified by this hook" in e.stdout.decode():
                all_lines = [line for line in e.stdout.decode().split("\n")]
                files = "\n".join(all_lines[3:])
                log.debug(f"The following files were modified by prettier:\n {files}")
            else:
                log.warning(
                    "There was an error running the prettier pre-commit hook.\n"
                    f"STDOUT: {e.stdout.decode()}\nSTDERR: {e.stderr.decode()}"
                )
    else:
        log.debug("Not in a git repository, skipping pre-commit hook.")


def dump_json_with_prettier(file_name, file_content):
    """Dump a JSON file and run prettier on it.
    Args:
        file_name (Path | str): A file identifier as a string or pathlib.Path.
        file_content (dict): Content to dump into the JSON file
    """
    with open(file_name, "w") as fh:
        json.dump(file_content, fh, indent=4)
    run_prettier_on_file(file_name)


def dump_yaml_with_prettier(file_name: Path | str, file_content: dict) -> None:
    """Dump a YAML file and run prettier on it.

    Args:
        file_name (Path | str): A file identifier as a string or pathlib.Path.
        file_content (dict): Content to dump into the YAML file
    """
    with open(file_name, "w") as fh:
        yaml.safe_dump(file_content, fh)
    run_prettier_on_file(file_name)


def ignore_file(lint_name: str, file_path: Path, dir_path: Path) -> list[list[str]]:
    """Ignore a file and add the result to the ignored list. Return the passed, failed, ignored and ignore_configs lists."""

    passed: list[str] = []
    failed: list[str] = []
    ignored: list[str] = []
    _, pipeline_conf = nf_core.utils.load_tools_config(dir_path)
    lint_conf = getattr(pipeline_conf, "lint", None) or None

    if lint_conf is None:
        ignore_entry: list[str] = []
    else:
        ignore_entry = lint_conf.get(lint_name, [])
    full_path = dir_path / file_path
    # Return a failed status if we can't find the file
    if not full_path.is_file():
        if isinstance(ignore_entry, bool) and not ignore_entry:
            ignored.append(f"`{file_path}` not found, but it is ignored.")
            ignore_entry = []
        else:
            failed.append(f"`{file_path}` not found.")
    else:
        passed.append(f"`{file_path}` found and not ignored.")

    # we handled the only case where ignore_entry should be a bool, convert it to a list, to make downstream code easier
    if isinstance(ignore_entry, bool):
        ignore_entry = []

    return [passed, failed, ignored, ignore_entry]
