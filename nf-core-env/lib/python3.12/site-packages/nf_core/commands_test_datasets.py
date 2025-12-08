import logging

import click
import rich

from nf_core.test_datasets.list import list_dataset_branches, list_datasets
from nf_core.test_datasets.search import search_datasets
from nf_core.utils import rich_force_colors

log = logging.getLogger(__name__)
stdout = rich.console.Console(force_terminal=rich_force_colors())


def test_datasets_list_branches(ctx: click.Context) -> None:
    """
    List all branches on the nf-core/test-datasets repository.
    Only lists test data and module test data based on the curated list
    of pipeline names [on the website](https://raw.githubusercontent.com/nf-core/website/refs/heads/main/public/pipeline_names.json).
    """
    list_dataset_branches()


def test_datasets_list_remote(ctx: click.Context, branch: str, generate_nf_path: bool, generate_dl_url: bool) -> None:
    """
    List all files on a given branch in the remote nf-core/testdatasets repository on github.
    The resulting files can be parsed as a nextflow path or a url for downloading.
    """
    list_datasets(branch, generate_nf_path, generate_dl_url)


def test_datasets_search(
    ctx: click.Context, branch: str, generate_nf_path: bool, generate_dl_url: bool, query: str | None
) -> None:
    """
    Search all files on a given branch in the remote nf-core/testdatasets repository on github
    with an interactive autocompleting prompt and print the file matching the query.
    Specifying a branch is required.
    The resulting file can optionally be parsed as a nextflow path or a url for downloading
    """
    search_datasets(branch, generate_nf_path=generate_nf_path, generate_dl_url=generate_dl_url, query=query)
