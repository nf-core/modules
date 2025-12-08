import json
import logging
from dataclasses import dataclass
from pathlib import Path

import questionary
import requests
import rich

from nf_core.utils import (
    determine_base_dir,
    fetch_wf_config,
    load_tools_config,
    nfcore_question_style,
    rich_force_colors,
)

log = logging.getLogger(__name__)
stdout = rich.console.Console(force_terminal=rich_force_colors())

# Name of nf-core/test-datasets github branch for modules
MODULES_BRANCH_NAME = "modules"


# Files / directories starting with one of the following in a git tree are ignored:
IGNORED_FILE_PREFIXES = [
    ".",
    "CITATION",
    "LICENSE",
    "README",
    "docs",
]


# Inline help text to display in autompletion prompts
AUTOCOMPLETION_HINT = "(press 'tab' to autocomplete)"


@dataclass
class GithubApiEndpoints:
    gh_api_base_url: str = "https://api.github.com"
    gh_orga: str = "nf-core"
    gh_repo: str = "test-datasets"

    def get_pipelines_list_url(self) -> str:
        return "https://raw.githubusercontent.com/nf-core/website/refs/heads/main/public/pipeline_names.json"

    def get_remote_tree_url_for_branch(self, branch: str) -> str:
        url = f"{self.gh_api_base_url}/repos/{self.gh_orga}/{self.gh_repo}/git/trees/{branch}?recursive=1"
        return url

    def get_file_download_url(self, branch: str, path: str) -> str:
        url = f"https://raw.githubusercontent.com/{self.gh_orga}/{self.gh_repo}/refs/heads/{branch}/{path}"
        return url


def get_remote_branch_names() -> list[str]:
    """
    List all branch names on the remote github repository for test-datasets for pipelines or modules.
    """
    try:
        url = GithubApiEndpoints().get_pipelines_list_url()
        response = requests.get(url)
        resp_json = response.json()

        if not response.ok:
            log.error(
                f"HTTP status code {response.status_code} received while fetching the list of branches at url: {response.url}"
            )
            return []

        branches = resp_json["pipeline"]  # pipeline branches from curated list
        branches += [MODULES_BRANCH_NAME]  # modules test-datasets branch

    except requests.exceptions.RequestException as e:
        log.error(f"Error while handling request to url {url}", e)
    except KeyError as e:
        log.error("Error parsing the list of branches received from Github API", e)
    except json.decoder.JSONDecodeError as e:
        log.error(f"Error parsing the list of branches received from Github API at url {response.url} as json", e)

    return branches


def get_remote_tree_for_branch(branch: str, only_files: bool = True, ignored_prefixes: list[str] = []) -> list[str]:
    """
    For a given branch name, return the file tree by querying the github API
    at the endpoint at `/repos/nf-core/test-datasets/git/trees/`
    """
    gh_filetree_file_value = "blob"  # value in nodes used to refer to "files"
    gh_response_filetree_key = "tree"  # key in response to refer to the filetree
    gh_filetree_type_key = "type"  # key in filetree nodes used to refer to their type
    gh_filetree_name_key = "path"  # key in filetree nodes used to refer to their name

    try:
        gh_api_url = GithubApiEndpoints(gh_repo="test-datasets")
        response = requests.get(gh_api_url.get_remote_tree_url_for_branch(branch))

        if not response.ok:
            log.error(
                f"HTTP status code {response.status_code} received while fetching the repository filetree at url {response.url}"
            )
            return []

        repo_tree = json.loads(response.text)[gh_response_filetree_key]

        if only_files:
            repo_tree = [node for node in repo_tree if node[gh_filetree_type_key] == gh_filetree_file_value]

        # filter by ignored_prefixes and extract names
        repo_files = []
        for node in repo_tree:
            for prefix in ignored_prefixes:
                if node[gh_filetree_name_key].startswith(prefix):
                    break
            else:
                repo_files.append(node[gh_filetree_name_key])

    except requests.exceptions.RequestException as e:
        log.error(f"Error while handling request to url {gh_api_url.get_remote_tree_url_for_branch(branch)}", e)

    except json.decoder.JSONDecodeError as e:
        log.error(f"Error parsing the repository filetree received from Github API at url {response.url} as json", e)

    return repo_files


def list_files_by_branch(
    branch: str = "",
    branches: list[str] = [],
    ignored_file_prefixes: list[str] = [
        ".",
        "CITATION",
        "LICENSE",
        "README",
        "docs",
    ],
) -> dict[str, list[str]]:
    """
    Lists files for all branches in the test-datasets github repo.
    Returns dictionary with branchnames as keys and file-lists as values
    """

    if len(branches) == 0:
        log.debug("Fetching list of remote branch names")
        branches = get_remote_branch_names()

    if branch:
        branches = list(filter(lambda b: b == branch, branches))
        if len(branches) == 0:
            log.error(f"No branches matching '{branch}'")

    log.debug("Fetching remote trees")
    tree = dict()
    for b in branches:
        tree[b] = get_remote_tree_for_branch(b, only_files=True, ignored_prefixes=ignored_file_prefixes)

    return tree


def create_pretty_nf_path(path: str, is_module_dataset: bool) -> str:
    """
    Generates a line of nexflow code with the full file path including a test data base path.
    """
    out = "params."
    out += "modules_" if is_module_dataset else "pipelines_"
    out += f'testdata_base_path + "{path}"'
    return out


def create_download_url(branch: str, path: str) -> str:
    """
    Generate a github API download url for the given path and branch
    """
    gh_api_url = GithubApiEndpoints(gh_repo="test-datasets")
    return gh_api_url.get_file_download_url(branch, path)


def get_or_prompt_branch(maybe_branch: str) -> tuple[str, list[str]]:
    """
    If branch is given, return a tuple of (maybe_branch, empty_list) else
    prompt the user to enter a branch name and return (branch_name, all_branches)
    """
    if maybe_branch:
        return (maybe_branch, [])

    else:
        all_branches = get_remote_branch_names()

        # Find pipeline / modules root directory
        base_dir: Path = determine_base_dir()

        # Read .nf-core.yml to identify repository_type
        _, tools_config = load_tools_config(base_dir)

        branch_prefill = ""
        # either modules or a pipeline branch
        if tools_config is not None:
            repo_type = tools_config.get("repository_type", None)
            if repo_type == MODULES_BRANCH_NAME:
                branch_prefill = MODULES_BRANCH_NAME
            elif repo_type == "pipeline":
                wf_config = fetch_wf_config(base_dir)
                pipeline_name = wf_config.get("manifest.name", "").split("/")[-1]
                if pipeline_name in all_branches:
                    branch_prefill = pipeline_name

        branch = questionary.autocomplete(
            "Branch name:",
            choices=sorted(all_branches),
            style=nfcore_question_style,
            default=branch_prefill,
            qmark=AUTOCOMPLETION_HINT,
        ).unsafe_ask()

        return branch, all_branches


def get_or_prompt_file_selection(files: list[str], query: str | None) -> str:
    """
    Prompt with autocompletion to enter a file from a list of files until a valid file is selected.
    """
    file_selected = False
    query = query if query is not None else ""  # ensure query is not None

    if query:
        # Check if only one file matches the query and directly return it
        filtered_files = [f for f in files if query in f]
        if len(filtered_files) == 1:
            selection = filtered_files[0]
            file_selected = True

    while not file_selected:
        selection = questionary.autocomplete(
            "File:", choices=files, style=nfcore_question_style, default=query, qmark=AUTOCOMPLETION_HINT
        ).unsafe_ask()

        file_selected = any([selection == file for file in files])
        if not file_selected:
            stdout.print("Please select a file.")

    return selection
