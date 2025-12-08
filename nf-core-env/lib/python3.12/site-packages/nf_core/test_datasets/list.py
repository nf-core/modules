import logging
import os

import rich

from nf_core.test_datasets.test_datasets_utils import (
    IGNORED_FILE_PREFIXES,
    MODULES_BRANCH_NAME,
    create_download_url,
    create_pretty_nf_path,
    get_or_prompt_branch,
    get_remote_branch_names,
    list_files_by_branch,
)
from nf_core.utils import rich_force_colors

stdout = rich.console.Console(force_terminal=rich_force_colors())
log = logging.getLogger(__name__)


def list_dataset_branches() -> None:
    """
    List all branches on the nf-core/test-datasets repository.
    Only lists test data and module test data based on the curated list
    of pipeline names [on the website](https://raw.githubusercontent.com/nf-core/website/refs/heads/main/public/pipeline_names.json).
    """
    remote_branches = get_remote_branch_names()

    table = rich.table.Table()
    table.add_column("Test-Dataset Branches")
    for b in remote_branches:
        table.add_row(b)
    stdout.print(table)


def list_datasets(
    maybe_branch: str = "",
    generate_nf_path: bool = False,
    generate_dl_url: bool = False,
    ignored_file_prefixes: list[str] = IGNORED_FILE_PREFIXES,
) -> None:
    """
    List all datasets for the given branch on the nf-core/test-datasets repository.
    If the given branch is empty or None, the user is prompted to enter one.

    Only lists test data and module test data based on the curated list
    of pipeline names [on the website](https://raw.githubusercontent.com/nf-core/website/refs/heads/main/public/pipeline_names.json).

    Ignores files with prefixes given in ignored_file_prefixes.
    """
    if generate_nf_path and generate_dl_url:
        log.warning("Ignoring url output as nextflow path output is enabled.")

    branch, all_branches = get_or_prompt_branch(maybe_branch)

    tree = list_files_by_branch(branch, all_branches, ignored_file_prefixes)

    out = []
    for b in tree.keys():
        files = sorted(tree[b])
        for f in files:
            if generate_nf_path:
                out.append(create_pretty_nf_path(f, branch == MODULES_BRANCH_NAME))
            elif generate_dl_url:
                out.append(create_download_url(branch, f))
            else:
                out.append(f)

    plain_text_output = generate_nf_path or generate_dl_url
    if plain_text_output:
        stdout.print(os.linesep.join(out))

    else:
        table = rich.table.Table()
        table.add_column("File", overflow="fold")

        for el in out:
            table.add_row(el)

        stdout.print(table)
