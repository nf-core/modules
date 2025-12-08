import logging

import rich.console
import rich.table

from nf_core.test_datasets.test_datasets_utils import (
    IGNORED_FILE_PREFIXES,
    MODULES_BRANCH_NAME,
    create_download_url,
    create_pretty_nf_path,
    get_or_prompt_branch,
    get_or_prompt_file_selection,
    list_files_by_branch,
)
from nf_core.utils import rich_force_colors

stdout = rich.console.Console(force_terminal=rich_force_colors())
log = logging.getLogger(__name__)


def search_datasets(
    maybe_branch: str = "",
    generate_nf_path: bool = False,
    generate_dl_url: bool = False,
    ignored_file_prefixes: list[str] = IGNORED_FILE_PREFIXES,
    plain_text_output: bool = False,
    query: str | None = "",
) -> None:
    """
    Search all files on a given branch in the remote nf-core/testdatasets repository on github
    with an interactive autocompleting prompt and print the file matching the query.

    Specifying a branch is required. If an empty string or None is specified as a branch,
    the user is prompted to selecte a branch as well.

    The resulting file can optionally be parsed as a nextflow path or a url for downloading
    """

    if generate_nf_path and generate_dl_url:
        log.warning("Ignoring url output as nextflow path output is enabled.")

    branch, all_branches = get_or_prompt_branch(maybe_branch)

    stdout.print("Searching files on branch: ", branch)
    tree = list_files_by_branch(branch, all_branches, ignored_file_prefixes)
    files = sum(tree.values(), [])  # flat representation of tree

    selection = get_or_prompt_file_selection(files, query)

    if generate_nf_path:
        stdout.print(create_pretty_nf_path(selection, branch == MODULES_BRANCH_NAME))
    elif generate_dl_url:
        stdout.print(create_download_url(branch, selection))
    elif plain_text_output:
        stdout.print(selection)
        stdout.print(create_pretty_nf_path(selection, branch == MODULES_BRANCH_NAME))
        stdout.print(create_download_url(branch, selection))
    else:
        table = rich.table.Table(show_header=False)
        table.add_column("")
        table.add_column("", overflow="fold")
        table.add_row("File Name:", selection)
        table.add_row("Nextflow Import:", create_pretty_nf_path(selection, branch == MODULES_BRANCH_NAME))
        table.add_row("Download Link:", create_download_url(branch, selection))
        stdout.print(table)
