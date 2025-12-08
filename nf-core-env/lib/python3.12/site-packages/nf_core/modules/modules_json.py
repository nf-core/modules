import copy
import datetime
import json
import logging
import os
import shutil
import tempfile
from pathlib import Path

import git
import questionary
import rich.prompt
from git.exc import GitCommandError
from typing_extensions import NotRequired, TypedDict  # for py<3.11

import nf_core.utils
from nf_core.components.components_utils import get_components_to_install
from nf_core.components.constants import NF_CORE_MODULES_NAME, NF_CORE_MODULES_REMOTE
from nf_core.modules.modules_repo import ModulesRepo
from nf_core.pipelines.lint_utils import dump_json_with_prettier

from ..components.components_differ import ComponentsDiffer

log = logging.getLogger(__name__)


class ModulesJsonModuleEntry(TypedDict):
    branch: str
    git_sha: str
    installed_by: list[str]
    patch: NotRequired[str]


class ModulesJsonType(TypedDict):
    name: str
    homePage: str
    repos: dict[str, dict[str, dict[str, dict[str, ModulesJsonModuleEntry]]]]


class ModulesJson:
    """
    An object for handling a 'modules.json' file in a pipeline
    """

    def __init__(self, pipeline_dir: str | Path) -> None:
        """
        Initialise the object.

        Args:
            pipeline_dir (str): The pipeline directory
        """
        self.directory = Path(pipeline_dir)
        self.modules_dir = self.directory / "modules"
        self.subworkflows_dir = self.directory / "subworkflows"
        self.modules_json_path = self.directory / "modules.json"
        self.modules_json: ModulesJsonType | None = None
        self.pipeline_modules = None
        self.pipeline_subworkflows = None
        self.pipeline_components: dict[str, list[tuple[str, str]]] | None = None

    def __str__(self):
        if self.modules_json is None:
            self.load()
        return json.dumps(self.modules_json, indent=4)

    def __repr__(self):
        return self.__str__()

    def create(self) -> None:
        """
        Creates the modules.json file from the modules and subworkflows installed in the pipeline directory

        Raises:
            UserWarning: If the creation fails
        """
        pipeline_config = nf_core.utils.fetch_wf_config(self.directory)
        pipeline_name = pipeline_config.get("manifest.name", "")
        pipeline_url = pipeline_config.get("manifest.homePage", "")
        new_modules_json = ModulesJsonType(name=pipeline_name, homePage=pipeline_url, repos={})

        if not self.modules_dir.exists():
            if rich.prompt.Confirm.ask(
                "[bold][blue]?[/] Can't find a ./modules directory. Would you like me to create one?", default=True
            ):
                log.info(f"Creating ./modules directory in '{self.directory}'")
                self.modules_dir.mkdir()
            else:
                raise UserWarning("Cannot proceed without a ./modules directory.")

        # Get repositories
        repos, _ = self.get_pipeline_module_repositories("modules", self.modules_dir)
        # Get all module/subworkflow names in the repos
        repo_module_names = self.get_component_names_from_repo(repos, self.modules_dir)
        repo_subworkflow_names = self.get_component_names_from_repo(repos, self.subworkflows_dir)

        # Add module/subworkflow info
        for repo_url, module_names, install_dir in sorted(repo_module_names):
            new_modules_json["repos"][repo_url] = {}
            new_modules_json["repos"][repo_url]["modules"] = {}
            new_modules_json["repos"][repo_url]["modules"][install_dir] = {}
            new_modules_json["repos"][repo_url]["modules"][install_dir] = self.determine_branches_and_shas(
                "modules", install_dir, repo_url, module_names
            )
        for repo_url, subworkflow_names, install_dir in sorted(repo_subworkflow_names):
            if repo_url not in new_modules_json["repos"]:  # Don't overwrite the repo if it was already added by modules
                new_modules_json["repos"][repo_url] = {}
            new_modules_json["repos"][repo_url]["subworkflows"] = {}
            new_modules_json["repos"][repo_url]["subworkflows"][install_dir] = {}
            new_modules_json["repos"][repo_url]["subworkflows"][install_dir] = self.determine_branches_and_shas(
                "subworkflows", install_dir, repo_url, subworkflow_names
            )

        # write the modules.json file and assign it to the object
        self.modules_json = new_modules_json
        self.dump()

    def get_component_names_from_repo(
        self, repos: dict[str, dict[str, dict[str, dict[str, dict[str, str | list[str]]]]]], directory: Path
    ) -> list[tuple[str, list[str], str]]:
        """
        Get component names from repositories in a pipeline.

        Args:
            repos (list): list of repository urls
            directory (str): modules directory or subworkflows directory

        Returns:
            [(str),[(str),(str)]]: list of tuples with repository url, component names and install directory
        """
        names = []
        for repo_url in repos:
            modules_repo = ModulesRepo(repo_url)
            if modules_repo is None:
                raise UserWarning(f"Could not find module repository for '{repo_url}' in '{directory}'")
            if modules_repo.repo_path is None:
                raise UserWarning(f"Could not find module repository path for '{repo_url}' in '{directory}'")
            components = (
                repo_url,
                [
                    str(Path(component_name).relative_to(directory / modules_repo.repo_path))
                    for component_name, _, file_names in os.walk(directory / modules_repo.repo_path)
                    if "main.nf" in file_names
                ],
                modules_repo.repo_path,
            )
            names.append(components)
        return names

    def get_pipeline_module_repositories(
        self, component_type: str, directory: Path, repos: dict | None = None
    ) -> tuple[dict[str, dict[str, dict[str, dict[str, dict[str, str | list[str]]]]]], dict[Path, Path]]:
        """
        Finds all module repositories in the modules and subworkflows directory.
        Ignores the local modules/subworkflows.

        Args:
            component_type (str): modules or subworkflows
            directory (Path): base directory for the module files
        Returns
            repos ([ (str, str, str) ]),
            renamed_dirs (dict[Path, Path]): List of tuples of repo name, repo
                                             remote URL and path to modules in
                                             repo
        """
        if repos is None:
            repos = {}
        # Check if there are any nf-core modules installed
        if (directory / NF_CORE_MODULES_NAME).exists() and NF_CORE_MODULES_REMOTE not in repos.keys():
            repos[NF_CORE_MODULES_REMOTE] = {}
        # The function might rename some directories, keep track of them
        renamed_dirs = {}
        # Check if there are any untracked repositories

        dirs_not_covered = self.dir_tree_uncovered(directory, [Path(ModulesRepo(url).repo_path) for url in repos])
        if len(dirs_not_covered) > 0:
            log.info(f"Found custom {component_type[:-1]} repositories when creating 'modules.json'")
            # Loop until all directories in the base directory are covered by a remote
            while len(dirs_not_covered) > 0:
                log.info(
                    "The following director{s} in the {t} directory are untracked: '{l}'".format(
                        s="ies" if len(dirs_not_covered) > 0 else "y",
                        t=component_type,
                        l="', '".join(str(dir.relative_to(directory)) for dir in dirs_not_covered),
                    )
                )
                nrepo_remote = questionary.text(
                    "Please provide a URL for for one of the repos contained in the untracked directories.",
                    style=nf_core.utils.nfcore_question_style,
                ).unsafe_ask()
                # Verify that the remote exists
                while True:
                    try:
                        git.Git().ls_remote(nrepo_remote)
                        break
                    except GitCommandError:
                        nrepo_remote = questionary.text(
                            "The provided remote does not seem to exist, please provide a new remote."
                        ).unsafe_ask()

                # Verify that there is a directory corresponding the remote
                nrepo_name = ModulesRepo(nrepo_remote).repo_path
                if nrepo_name is None:
                    raise UserWarning(f"Could not find the repository name for '{nrepo_remote}'")
                if not (directory / nrepo_name).exists():
                    log.info(
                        "The provided remote does not seem to correspond to a local directory. "
                        "The directory structure should be the same as in the remote."
                    )
                    dir_name = questionary.text(
                        "Please provide the correct directory, it will be renamed. If left empty, the remote will be ignored.",
                        style=nf_core.utils.nfcore_question_style,
                    ).unsafe_ask()
                    if dir_name:
                        old_path = directory / dir_name
                        new_path = directory / nrepo_name
                        old_path.rename(new_path)
                        renamed_dirs[old_path] = new_path
                    else:
                        continue

                if nrepo_remote not in repos:
                    repos[nrepo_remote] = {}
                if component_type not in repos[nrepo_remote]:
                    repos[nrepo_remote][component_type] = {}
                repos[nrepo_remote][component_type][nrepo_name] = {}
                dirs_not_covered = self.dir_tree_uncovered(
                    directory, [Path(name) for url in repos for name in repos[url][component_type]]
                )

        return repos, renamed_dirs

    def dir_tree_uncovered(self, components_directory, repos):
        """
        Does a BFS of the modules/subworkflos directory to look for directories that
        are not tracked by a remote. The 'repos' argument contains the
        directories that are currently covered by remote, and it and its
        subdirectories are therefore ignore.

        Args:
            components_directory (Path): Base path of modules or subworkflows in pipeline
            repos ([ Path ]): List of repos that are covered by a remote

        Returns:
            dirs_not_covered ([ Path ]): A list of directories that are currently not covered by any remote.
        """
        # Initialise the FIFO queue. Note that we assume the directory to be correctly
        # configured, i.e. no files etc.
        fifo = [subdir for subdir in components_directory.iterdir() if subdir.stem != "local"]
        depth = 1
        dirs_not_covered = []
        while len(fifo) > 0:
            temp_queue = []
            repos_at_level = {Path(*repo.parts[:depth]): len(repo.parts) for repo in repos}
            for directory in fifo:
                rel_dir = directory.relative_to(components_directory)
                if rel_dir in repos_at_level.keys():
                    # Go the next depth if this directory is not one of the repos
                    if depth < repos_at_level[rel_dir]:
                        temp_queue.extend(directory.iterdir())
                else:
                    # Otherwise add the directory to the ones not covered
                    dirs_not_covered.append(directory)
            fifo = temp_queue
            depth += 1
        return dirs_not_covered

    def determine_branches_and_shas(
        self,
        component_type: str,
        install_dir: str | Path,
        remote_url: str,
        components: list[str],
    ) -> dict[str, ModulesJsonModuleEntry]:
        """
        Determines what branch and commit sha each module/subworkflow in the pipeline belongs to

        Assumes all modules/subworkflows are installed from the default branch. If it fails to find the
        module/subworkflow in the default branch, it prompts the user with the available branches

        Args:
            install_dir (str): The name of the directory inside modules or subworkflows where components are installed
            remote_url (str): The url to the remote repository
            components ([str]): List of names of installed modules/subworkflows from the repository

        Returns:
            (dict[str, dict[str, str]]): The module.json entries for the modules/subworkflows
                                         from the repository
        """
        default_modules_repo = ModulesRepo(remote_url=remote_url)
        if component_type == "modules":
            repo_path = self.modules_dir / install_dir
        elif component_type == "subworkflows":
            repo_path = self.subworkflows_dir / install_dir
        else:
            raise ValueError(f"Unknown component type '{component_type}'")
        # Get the branches present in the repository, as well as the default branch
        available_branches = ModulesRepo.get_remote_branches(remote_url)
        sb_local = []
        dead_components = []
        repo_entry: dict[str, ModulesJsonModuleEntry] = {}
        for component in sorted(components):
            modules_repo = default_modules_repo
            component_path = Path(repo_path, component)
            correct_commit_sha = None
            tried_branches = {default_modules_repo.branch}
            found_sha = False
            while True:
                # If the module/subworkflow is patched
                patch_file = component_path / f"{component}.diff"
                if patch_file.is_file():
                    temp_module_dir = self.try_apply_patch_reverse(
                        component_type, component, install_dir, patch_file, component_path
                    )
                    correct_commit_sha = self.find_correct_commit_sha(
                        component_type, component, temp_module_dir, modules_repo
                    )
                else:
                    correct_commit_sha = self.find_correct_commit_sha(
                        component_type, component, component_path, modules_repo
                    )
                    if correct_commit_sha is None:
                        # Check in the old path
                        correct_commit_sha = self.find_correct_commit_sha(
                            component_type, component, repo_path / component_type / component, modules_repo
                        )
                if correct_commit_sha is None:
                    log.info(
                        f"Was unable to find matching {component_type[:-1]} files in the {modules_repo.branch} branch."
                    )
                    choices = [{"name": "No", "value": False}] + [
                        {"name": branch, "value": branch} for branch in (available_branches - tried_branches)
                    ]
                    branch = questionary.select(
                        f"Was the {component_type[:-1]} '{component}' installed from a different branch in the remote?\nSelect 'No' for a local {component_type[:-1]}",
                        choices=choices,
                        style=nf_core.utils.nfcore_question_style,
                    ).unsafe_ask()
                    if not branch:
                        action = questionary.select(
                            f"{component_type[:-1].title()} is untracked '{component}'. Please select what action to take",
                            choices=[
                                {"name": "Move the directory to 'local'", "value": 0},
                                {"name": "Remove the files", "value": 1},
                            ],
                            style=nf_core.utils.nfcore_question_style,
                        ).unsafe_ask()
                        if action == 0:
                            sb_local.append(component)
                        else:
                            dead_components.append(component)
                        break
                    # Create a new modules repo with the selected branch, and retry find the sha
                    modules_repo = ModulesRepo(remote_url=remote_url, branch=branch, no_pull=True, hide_progress=True)
                else:
                    found_sha = True
                    break
            if found_sha and correct_commit_sha is not None:
                repo_entry[component] = {
                    "branch": modules_repo.branch,
                    "git_sha": correct_commit_sha,
                    "installed_by": [component_type],
                }

        # Clean up the modules/subworkflows we were unable to find the sha for
        for component in sb_local:
            log.debug(f"Moving {component_type[:-1]} '{Path(install_dir, component)}' to 'local' directory")
            self.move_component_to_local(component_type, component, str(install_dir))

        for component in dead_components:
            log.debug(f"Removing {component_type[:-1]} {Path(install_dir, component)}'")
            shutil.rmtree(repo_path / component)

        return repo_entry

    def find_correct_commit_sha(
        self,
        component_type: str,
        component_name: str | Path,
        component_path: str | Path,
        modules_repo: ModulesRepo,
    ) -> str | None:
        """
        Returns the SHA for the latest commit where the local files are identical to the remote files
        Args:
            component_type (str): modules or subworkflows
            component_name (str): Name of module/subowrkflow
            component_path (str): Path to module/subworkflow in local repo
            modules_repo (str): Remote repo for module/subworkflow
        Returns:
            commit_sha (str): The latest commit SHA where local files are identical to remote files,
                              or None if no commit is found
        """
        # Find the correct commit SHA for the local module/subworkflow files.
        # We iterate over the commit history for the module/subworkflow until we find
        # a revision that matches the file contents
        commit_shas = (
            commit["git_sha"]
            for commit in modules_repo.get_component_git_log(component_name, component_type, depth=1000)
        )
        for commit_sha in commit_shas:
            if all(
                modules_repo.component_files_identical(
                    component_name, component_path, commit_sha, component_type
                ).values()
            ):
                return commit_sha
        return None

    def move_component_to_local(self, component_type: str, component: str, repo_name: str):
        """
        Move a module/subworkflow to the 'local' directory

        Args:
            component_type (str): The type of component, either 'modules' or 'subworkflows'
            component (str): The name of the module/subworkflow
            repo_name (str): The name of the repository the module resides in
        """
        if component_type == "modules":
            directory = self.modules_dir
        elif component_type == "subworkflows":
            directory = self.subworkflows_dir
        else:
            raise ValueError(f"Unknown component type '{component_type}'")
        current_path = directory / repo_name / component
        local_dir = directory / "local"
        if not local_dir.exists():
            local_dir.mkdir(parents=True)

        to_name = component
        # Check if there is already a subdirectory with the name
        while (local_dir / to_name).exists():
            # Add a time suffix to the path to make it unique
            # (do it again and again if it didn't work out...)
            to_name += f"-{datetime.datetime.now().strftime('%y%m%d%H%M%S')}"
        shutil.move(str(current_path), local_dir / to_name)

    def unsynced_components(self) -> tuple[list[str], list[str], dict]:
        """
        Compute the difference between the modules/subworkflows in the directory and the
        modules/subworkflows in the 'modules.json' file. This is done by looking at all
        directories containing a 'main.nf' file

        Returns:
            (untrack_dirs ([ str ]), missing_installation (dict)): Directories that are not tracked
            by the modules.json file, and modules/subworkflows in the modules.json where
            the installation directory is missing
        """
        assert self.modules_json is not None  # mypy
        # Add all modules from modules.json to missing_installation
        missing_installation = copy.deepcopy(self.modules_json["repos"])
        # Obtain the path of all installed modules
        module_dirs = [
            Path(dir_name).relative_to(self.modules_dir)
            for dir_name, _, file_names in os.walk(self.modules_dir)
            if "main.nf" in file_names and not str(Path(dir_name).relative_to(self.modules_dir)).startswith("local")
        ]
        untracked_dirs_modules, missing_installation = self.parse_dirs(module_dirs, missing_installation, "modules")

        # Obtain the path of all installed subworkflows
        subworkflow_dirs = [
            Path(dir_name).relative_to(self.subworkflows_dir)
            for dir_name, _, file_names in os.walk(self.subworkflows_dir)
            if "main.nf" in file_names
            and not str(Path(dir_name).relative_to(self.subworkflows_dir)).startswith("local")
        ]
        untracked_dirs_subworkflows, missing_installation = self.parse_dirs(
            subworkflow_dirs, missing_installation, "subworkflows"
        )

        return untracked_dirs_modules, untracked_dirs_subworkflows, missing_installation

    def parse_dirs(self, dirs: list[Path], missing_installation: dict, component_type: str) -> tuple[list[str], dict]:
        """
        Parse directories and check if they are tracked in the modules.json file

        Args:
            dirs ([ Path ]): List of directories to check
            missing_installation (dict): Dictionary with the modules.json entries
            component_type (str): The type of component, either 'modules' or 'subworkflows'

        Returns:
            (untracked_dirs ([ Path ]), missing_installation (dict)): List of directories that are not tracked
            by the modules.json file, and the updated missing_installation dictionary
        """

        untracked_dirs = []
        for dir_ in dirs:
            # Check if the module/subworkflows directory exists in modules.json
            install_dir = dir_.parts[0]
            component = "/".join(dir_.parts[1:])
            component_in_file = False
            git_url = ""

            for repo in missing_installation:
                if component_type in missing_installation[repo]:
                    if install_dir in missing_installation[repo][component_type]:
                        if component in missing_installation[repo][component_type][install_dir]:
                            component_in_file = True
                            git_url = repo
                            break
            if not component_in_file:
                # If it is not, add it to the list of missing components
                untracked_dirs.append(component)
            else:
                # If it does, remove the component from missing_installation
                module_repo = missing_installation[git_url]
                # Check if the entry has a git sha and branch before removing
                components_dict = module_repo[component_type][install_dir]
                if "git_sha" not in components_dict[component] or "branch" not in components_dict[component]:
                    self.determine_branches_and_shas(component_type, component, git_url, [component])
                # Remove the module/subworkflow from modules/subworkflows without installation
                module_repo[component_type][install_dir].pop(component)
                if len(module_repo[component_type][install_dir]) == 0:
                    # If no modules/subworkflows with missing installation left, remove the install_dir from missing_installation
                    missing_installation[git_url][component_type].pop(install_dir)
                if len(module_repo[component_type]) == 0:
                    # If no modules/subworkflows with missing installation left, remove the component_type from missing_installation
                    missing_installation[git_url].pop(component_type)
                if len(module_repo) == 0:
                    # If no modules/subworkflows with missing installation left, remove the git_url from missing_installation
                    missing_installation.pop(git_url)

        return untracked_dirs, missing_installation

    def has_git_url_and_modules(self) -> bool:
        """
        Check that all repo entries in the modules.json
        has a git url and a modules dict entry
        Returns:
            (bool): True if they are found for all repos, False otherwise
        """
        assert self.modules_json is not None  # mypy
        for repo_url, repo_entry in self.modules_json.get("repos", {}).items():
            if "modules" not in repo_entry:
                if "subworkflows" in repo_entry:
                    continue
                log.warning(f"modules.json entry {repo_entry} does not have a modules entry")
                return False
            elif (
                not isinstance(repo_url, str)
                or repo_url == ""
                or not (
                    repo_url.startswith("http")
                    or repo_url.startswith("ftp")
                    or repo_url.startswith("ssh")
                    or repo_url.startswith("git")
                )
                or not isinstance(repo_entry["modules"], dict)
                or repo_entry["modules"] == {}
            ):
                log.debug(f"modules.json entry {repo_entry} has non-string or empty entries for git_url or modules.")
                return False
        return True

    def reinstall_repo(self, install_dir, remote_url, module_entries):
        """
        Reinstall modules from a repository

        Args:
            install_dir (str): The name of directory where modules are installed
            remote_url (str): The git url of the remote repository
            module_entries ([ dict[str, dict[str, str]] ]): Module entries with
            branch and git sha info

        Returns:
            ([ str ]): List of modules that we failed to install
        """
        branches_and_mods = {}
        failed_to_install = []
        for module, module_entry in module_entries.items():
            if "git_sha" not in module_entry or "branch" not in module_entry:
                failed_to_install.append(module)
            else:
                branch = module_entry["branch"]
                sha = module_entry["git_sha"]
                if branch not in branches_and_mods:
                    branches_and_mods[branch] = []
                branches_and_mods[branch].append((module, sha))

        for branch, modules in branches_and_mods.items():
            try:
                modules_repo = ModulesRepo(remote_url=remote_url, branch=branch)
            except LookupError as e:
                log.error(e)
                failed_to_install.extend(modules)
            for module, sha in modules:
                if not modules_repo.install_component(module, self.modules_dir / install_dir, sha, "modules"):
                    log.warning(
                        f"Could not install module '{Path(self.modules_dir, install_dir, module)}' - removing from modules.json"
                    )
                    failed_to_install.append(module)
        return failed_to_install

    def check_up_to_date(self):
        """
        Checks whether the modules and subworkflows installed in the directory
        are consistent with the entries in the 'modules.json' file and vice versa.

        If a module/subworkflow has an entry in the 'modules.json' file but is missing in the directory,
        we first try to reinstall the module/subworkflow from the remote and if that fails we remove the entry
        in 'modules.json'.

        If a module/subworkflow is installed but the entry in 'modules.json' is missing we iterate through
        the commit log in the remote to try to determine the SHA.

        Check that we have the "installed_by" value in 'modules.json', otherwise add it.
        Assume that the modules/subworkflows were installed by an nf-core command (don't track installed by subworkflows).
        """
        dump_modules_json = False
        try:
            self.load()
            if not self.has_git_url_and_modules():
                raise UserWarning

            assert self.modules_json is not None  # mypy
            # check that all "installed_by" entries are lists and not strings
            # [these strings come from an older dev version, so this check can probably be removed in a future release]
            for _, repo_entry in self.modules_json.get("repos", {}).items():
                for component_type in ["modules", "subworkflows"]:
                    if component_type in repo_entry:
                        for install_dir, install_dir_entry in repo_entry[component_type].items():
                            for _, component in install_dir_entry.items():
                                if "installed_by" in component and isinstance(component["installed_by"], str):
                                    log.debug(f"Updating {component} in modules.json")
                                    dump_modules_json = True
                                    component["installed_by"] = [component["installed_by"]]
        except UserWarning:
            log.info("The 'modules.json' file is not up to date. Recreating the 'modules.json' file.")
            dump_modules_json = True
            self.create()

        # Get unsynced components
        (
            modules_missing_from_modules_json,
            subworkflows_missing_from_modules_json,
            missing_installation,
        ) = self.unsynced_components()

        # If there are any modules/subworkflows left in 'modules.json' after all installed are removed,
        # we try to reinstall them
        if len(missing_installation) > 0:
            if "subworkflows" in [
                c_type for _, repo_content in missing_installation.items() for c_type in repo_content.keys()
            ]:
                self.resolve_missing_installation(missing_installation, "subworkflows")
            if "modules" in [
                c_type for _, repo_content in missing_installation.items() for c_type in repo_content.keys()
            ]:
                self.resolve_missing_installation(missing_installation, "modules")

        # If some modules/subworkflows didn't have an entry in the 'modules.json' file
        # we try to determine the SHA from the commit log of the remote
        if len(modules_missing_from_modules_json) > 0:
            dump_modules_json = True
            self.resolve_missing_from_modules_json(modules_missing_from_modules_json, "modules")
        if len(subworkflows_missing_from_modules_json) > 0:
            dump_modules_json = True
            self.resolve_missing_from_modules_json(subworkflows_missing_from_modules_json, "subworkflows")
        assert self.modules_json is not None  # mypy
        # If the "installed_by" value is not present for modules/subworkflows, add it.
        for repo, repo_content in self.modules_json["repos"].items():
            for component_type, dir_content in repo_content.items():
                for install_dir, installed_components in dir_content.items():
                    for component, component_features in installed_components.items():
                        if "installed_by" not in component_features:
                            dump_modules_json = True
                            self.modules_json["repos"][repo][component_type][install_dir][component]["installed_by"] = [
                                component_type
                            ]

        # Recreate "installed_by" entry
        original_pipeline_components = self.pipeline_components
        self.pipeline_components = None
        subworkflows_dict = self.get_all_components("subworkflows")
        if subworkflows_dict:
            dump_modules_json = True
            for repo, subworkflows in subworkflows_dict.items():
                for org, subworkflow in subworkflows:
                    self.recreate_dependencies(repo, org, {"name": subworkflow})
        self.pipeline_components = original_pipeline_components

        if dump_modules_json:
            self.dump(run_prettier=True)
        return True

    def load(self) -> None:
        """
        Loads the modules.json file into the variable 'modules_json'

        Sets the modules_json attribute to the loaded file.

        Raises:
            UserWarning: If the modules.json file is not found
        """
        try:
            with open(self.modules_json_path) as fh:
                try:
                    self.modules_json = json.load(fh)
                except json.JSONDecodeError as e:
                    raise UserWarning(f"Unable to load JSON file '{self.modules_json_path}' due to error {e}")

        except FileNotFoundError:
            raise UserWarning("File 'modules.json' is missing")

    def update(
        self,
        component_type: str,
        modules_repo: ModulesRepo,
        component_name: str,
        component_version: str,
        installed_by: list[str] | None,
        installed_by_log: list[str] | None = None,
        write_file: bool = True,
    ) -> bool:
        """
        Updates the 'module.json' file with new module/subworkflow info

        Args:
            component_type (str): modules or subworkflows
            modules_repo (ModulesRepo): A ModulesRepo object configured for the new module/subworkflow
            component_name (str): Name of new module/subworkflow
            component_version (str): git SHA for the new module/subworkflow entry
            installed_by_log (list): previous tracing of installed_by that needs to be added to 'modules.json'
            write_file (bool): whether to write the updated modules.json to a file.

        Returns:
            bool: True if the module/subworkflow was successfully added to the 'modules.json' file
        """
        if installed_by_log is None:
            installed_by_log = []

        if self.modules_json is None:
            self.load()
            assert self.modules_json is not None  # mypy
        repo_name = modules_repo.repo_path
        remote_url = modules_repo.remote_url
        branch = modules_repo.branch

        if remote_url not in self.modules_json["repos"]:
            self.modules_json["repos"][remote_url] = {component_type: {repo_name: {}}}
        if component_type not in self.modules_json["repos"][remote_url]:
            self.modules_json["repos"][remote_url][component_type] = {repo_name: {}}
        repo_component_entry = self.modules_json["repos"][remote_url][component_type][repo_name]
        if component_name not in repo_component_entry:
            repo_component_entry[component_name] = {"branch": "", "git_sha": "", "installed_by": []}
        repo_component_entry[component_name]["git_sha"] = component_version
        repo_component_entry[component_name]["branch"] = branch
        try:
            if installed_by not in repo_component_entry[component_name]["installed_by"] and installed_by is not None:
                repo_component_entry[component_name]["installed_by"] += installed_by
        finally:
            new_installed_by = repo_component_entry[component_name]["installed_by"] + list(installed_by_log)
            repo_component_entry[component_name]["installed_by"] = sorted([*set(new_installed_by)])

        # Sort the 'modules.json' repo entries
        self.modules_json["repos"] = nf_core.utils.sort_dictionary(self.modules_json["repos"])
        if write_file:
            self.dump()
        return True

    def remove_entry(self, component_type, name, repo_url, install_dir, removed_by=None):
        """
        Removes an entry from the 'modules.json' file.

        Args:
            component_type (str): Type of component [modules, subworkflows]
            name (str): Name of the component to be removed
            repo_url (str): URL of the repository containing the component
            install_dir (str): Name of the directory where components are installed
            removed_by (str): Name of the component that wants to remove the component
        Returns:
            (bool): return True if the component was removed, False if it was not found or is still depended on
        """

        if removed_by is None or removed_by == name:
            removed_by = component_type
        if not self.modules_json:
            return False
        if repo_url in self.modules_json.get("repos", {}):
            repo_entry = self.modules_json["repos"][repo_url]
            if name in repo_entry[component_type].get(install_dir, {}):
                if removed_by in repo_entry[component_type][install_dir][name]["installed_by"]:
                    self.modules_json["repos"][repo_url][component_type][install_dir][name]["installed_by"].remove(
                        removed_by
                    )
                    # clean up empty entries
                    if len(repo_entry[component_type][install_dir][name]["installed_by"]) == 0:
                        self.modules_json["repos"][repo_url][component_type][install_dir].pop(name)
                        if len(repo_entry[component_type][install_dir]) == 0:
                            self.modules_json["repos"][repo_url].pop(component_type)
                            if len(repo_entry) == 0:
                                self.modules_json["repos"].pop(repo_url)
                        # write the updated modules.json file
                        self.dump()
                        return True
                    self.dump()
                    return False
            else:
                log.warning(
                    f"{component_type[:-1].title()} '{install_dir}/{name}' is missing from 'modules.json' file."
                )
                return False

        else:
            log.warning(f"{component_type[:-1].title()} '{install_dir}/{name}' is missing from 'modules.json' file.")
            return False

        return False

    def add_patch_entry(self, component_type, component_name, repo_url, install_dir, patch_filename, write_file=True):
        """
        Adds (or replaces) the patch entry for a module
        """
        if self.modules_json is None:
            self.load()
            assert self.modules_json is not None  # mypy

        if repo_url not in self.modules_json["repos"]:
            raise LookupError(f"Repo '{repo_url}' not present in 'modules.json'")
        if component_name not in self.modules_json["repos"][repo_url][component_type][install_dir]:
            raise LookupError(
                f"{component_type[:-1].title()} '{install_dir}/{component_name}' not present in 'modules.json'"
            )
        self.modules_json["repos"][repo_url][component_type][install_dir][component_name]["patch"] = str(patch_filename)
        if write_file:
            self.dump()

    def remove_patch_entry(self, module_name, repo_url, install_dir, write_file=True):
        if self.modules_json is None:
            self.load()
            assert self.modules_json is not None  # mypy

        try:
            del self.modules_json["repos"][repo_url]["modules"][install_dir][module_name]["patch"]
        except KeyError:
            log.warning("No patch entry in 'modules.json' to remove")
        if write_file:
            self.dump()

    def get_patch_fn(self, component_type, component_name, repo_url, install_dir):
        """
        Get the patch filename of a component

        Args:
            component_name (str): The name of the component
            repo_url (str): The URL of the repository containing the component
            install_dir (str): The name of the directory where components are installed

        Returns:
            (str): The patch filename for the component, None if not present
        """
        if self.modules_json is None:
            self.load()
            assert self.modules_json is not None  # mypy
        path = (
            self.modules_json["repos"]
            .get(repo_url, {})
            .get(component_type)
            .get(install_dir)
            .get(component_name, {})
            .get("patch")
        )
        return Path(path) if path is not None else None

    def try_apply_patch_reverse(self, component_type, component, repo_name, patch_relpath, component_dir):
        """
        Try reverse applying a patch file to the modified module or subworkflow files

        Args:
            component_type (str): The type of component [modules, subworkflows]
            component (str): The name of the module or subworkflow
            repo_name (str): The name of the repository where the component resides
            patch_relpath (Path | str): The path to patch file in the pipeline
            component_dir (Path | str): The component directory in the pipeline

        Returns:
            (Path | str): The path of the folder where the component patched files are

        Raises:
            LookupError: If patch was not applied
        """
        component_fullname = str(Path(repo_name, component))
        patch_path = Path(self.directory / patch_relpath)

        try:
            new_files = ComponentsDiffer.try_apply_patch(
                component_type, component, repo_name, patch_path, component_dir, reverse=True
            )
        except LookupError as e:
            raise LookupError(
                f"Failed to apply patch in reverse for {component_type[:-1]} '{component_fullname}' due to: {e}"
            )

        # Write the patched files to a temporary directory
        log.debug("Writing patched files to tmpdir")
        temp_dir = Path(tempfile.mkdtemp())
        temp_component_dir = temp_dir / component
        temp_component_dir.mkdir(parents=True, exist_ok=True)
        for file, new_content in new_files.items():
            fn = temp_component_dir / file
            with open(fn, "w") as fh:
                fh.writelines(new_content)

        return temp_component_dir

    def repo_present(self, repo_name):
        """
        Checks if a repo is present in the modules.json file
        Args:
            repo_name (str): Name of the repository
        Returns:
            (bool): Whether the repo exists in the modules.json
        """
        if self.modules_json is None:
            self.load()
            assert self.modules_json is not None  # mypy

        return repo_name in self.modules_json.get("repos", {})

    def component_present(self, module_name, repo_url, install_dir, component_type):
        """
        Checks if a module is present in the modules.json file
        Args:
            module_name (str): Name of the module
            repo_url (str): URL of the repository
            install_dir (str): Name of the directory where modules are installed
            component_type (str): Type of component [modules, subworkflows]
        Returns:
            (bool): Whether the module is present in the 'modules.json' file
        """
        if self.modules_json is None:
            self.load()
            assert self.modules_json is not None  # mypy
        return module_name in self.modules_json.get("repos", {}).get(repo_url, {}).get(component_type, {}).get(
            install_dir, {}
        )

    def get_modules_json(self) -> ModulesJsonType:
        """
        Returns a copy of the loaded modules.json

        Returns:
            (dict): A copy of the loaded modules.json
        """
        if self.modules_json is None:
            self.load()
            assert self.modules_json is not None  # mypy
        return copy.deepcopy(self.modules_json)

    def get_component_version(self, component_type, component_name, repo_url, install_dir):
        """
        Returns the version of a module or subworkflow

        Args:
            component_name (str): Name of the module/subworkflow
            repo_url (str): URL of the repository
            install_dir (str): Name of the directory where modules/subworkflows are installed

        Returns:
            (str): The git SHA of the module/subworkflow if it exists, None otherwise
        """
        if self.modules_json is None:
            self.load()
            assert self.modules_json is not None  # mypy
        return (
            self.modules_json.get("repos", {})
            .get(repo_url, {})
            .get(component_type, {})
            .get(install_dir, {})
            .get(component_name, {})
            .get("git_sha", None)
        )

    def get_module_version(self, module_name: str, repo_url: str, install_dir: str) -> str | None:
        """
        Returns the version of a module

        Args:
            module_name (str): Name of the module
            repo_url (str): URL of the repository
            install_dir (str): Name of the directory where modules are installed

        Returns:
            (str): The git SHA of the module if it exists, None otherwise
        """
        if self.modules_json is None:
            self.load()
            assert self.modules_json is not None  # mypy
        try:
            sha = self.modules_json["repos"][repo_url]["modules"][install_dir][module_name]["git_sha"]
        except KeyError:
            sha = None
        return sha

    def get_subworkflow_version(self, subworkflow_name, repo_url, install_dir):
        """
        Returns the version of a subworkflow

        Args:
            subworkflow_name (str): Name of the module
            repo_url (str): URL of the repository
            install_dir (str): Name of the directory where subworkflows are installed

        Returns:
            (str): The git SHA of the subworkflow if it exists, None otherwise
        """
        if self.modules_json is None:
            self.load()
            assert self.modules_json is not None  # mypy
        return (
            self.modules_json.get("repos", {})
            .get(repo_url, {})
            .get("subworkflows", {})
            .get(install_dir, {})
            .get(subworkflow_name, {})
            .get("git_sha", None)
        )

    def get_all_components(self, component_type: str) -> dict[str, list[tuple[(str, str)]]]:
        """
        Retrieves all pipeline modules/subworkflows that are reported in the modules.json

        Returns:
            (dict[str, [(str, str)]]): Dictionary indexed with the repo urls, with a
                                list of tuples (component_dir, components) as values
        """
        if self.modules_json is None:
            self.load()
            assert self.modules_json is not None  # mypy

        if self.pipeline_components is None:
            self.pipeline_components = {}
            for repo, repo_entry in self.modules_json.get("repos", {}).items():
                if component_type in repo_entry:
                    for directory, components in repo_entry[component_type].items():
                        self.pipeline_components[repo] = [(directory, m) for m in components]

        return self.pipeline_components

    def get_dependent_components(
        self,
        component_type,
        name,
        dependent_components,
    ) -> dict[str, tuple[str, str, str]]:
        """
        Retrieves all pipeline modules/subworkflows that are reported in the modules.json
        as being installed by the given component

        Args:
            component_type (str): Type of component [modules, subworkflows]
            name (str): Name of the component to find dependencies for

        Returns:
            (dict[str: str,]): Dictionary indexed with the component names, with component_type as value
        """

        if self.modules_json is None:
            self.load()
            assert self.modules_json is not None  # mypy
        component_types = ["modules"] if component_type == "modules" else ["modules", "subworkflows"]
        # Find all components that have an entry of install by of  a given component, recursively call this function for subworkflows
        for type in component_types:
            for repo_url in self.modules_json["repos"].keys():
                modules_repo = ModulesRepo(repo_url)
                install_dir = modules_repo.repo_path
                try:
                    for comp in self.modules_json["repos"][repo_url][type][install_dir]:
                        if name in self.modules_json["repos"][repo_url][type][install_dir][comp]["installed_by"]:
                            dependent_components[comp] = (repo_url, install_dir, type)
                except KeyError as e:
                    # This exception will raise when there are only modules installed
                    log.debug(f"Trying to retrieve all {type}. There aren't {type} installed. Failed with error {e}")
                    continue

        return dependent_components

    def get_installed_by_entries(self, component_type, name):
        """
        Retrieves all entries of installed_by for a given component

        Args:
            component_type (str): Type of component [modules, subworkflows]
            name (str): Name of the component to find dependencies for

        Returns:
            (list): The list of installed_by entries

        """
        if self.modules_json is None:
            self.load()
            assert self.modules_json is not None  # mypy
        installed_by_entries = {}
        for _, repo_entry in self.modules_json.get("repos", {}).items():
            if component_type in repo_entry:
                for _, components in repo_entry[component_type].items():
                    if name in components:
                        installed_by_entries = components[name]["installed_by"]
                        break

        return installed_by_entries

    def get_component_branch(self, component_type: str, component: str, repo_url: str, install_dir: str) -> str:
        """
        Gets the branch from which the module/subworkflow was installed

        Returns:
            (str): The branch name
        Raises:
            LookupError: If there is no branch entry in the `modules.json`
        """
        if self.modules_json is None:
            self.load()
            assert self.modules_json is not None  # mypy
        try:
            branch = self.modules_json["repos"][repo_url][component_type][install_dir][component]["branch"]
        except (KeyError, TypeError):
            branch = None
        if branch is None:
            raise LookupError(
                f"Could not find branch information for component '{Path(install_dir, component)}'."
                f"Please remove the 'modules.json' and rerun the command to recreate it"
            )
        return branch

    def dump(self, run_prettier: bool = False) -> None:
        """
        Sort the modules.json, and write it to file
        """
        # Sort the modules.json
        if self.modules_json is None:
            self.load()
        if self.modules_json is not None:
            self.modules_json["repos"] = nf_core.utils.sort_dictionary(self.modules_json["repos"])
            if run_prettier:
                dump_json_with_prettier(self.modules_json_path, self.modules_json)
            else:
                with open(self.modules_json_path, "w") as fh:
                    json.dump(self.modules_json, fh, indent=4)

    def resolve_missing_installation(self, missing_installation: dict, component_type: str) -> None:
        missing_but_in_mod_json = [
            f"'{component_type}/{install_dir}/{component}'"
            for repo_url, contents in missing_installation.items()
            for install_dir, dir_contents in contents[component_type].items()
            for component in dir_contents
        ]
        log.info(
            f"Reinstalling {component_type} found in 'modules.json' but missing from directory: {', '.join(missing_but_in_mod_json)}"
        )

        remove_from_mod_json = {}
        for repo_url, contents in missing_installation.items():
            for install_dir, component_entries in contents[component_type].items():
                remove_from_mod_json[(repo_url, install_dir)] = self.reinstall_repo(
                    install_dir, repo_url, component_entries
                )

        # If the reinstall fails, we remove those entries in 'modules.json'
        if sum(map(len, remove_from_mod_json.values())) > 0:
            uninstallable_components = [
                f"'{install_dir}/{component}'"
                for (repo_url, install_dir), components in remove_from_mod_json.items()
                for component in components
            ]
            if len(uninstallable_components) == 1:
                log.info(f"Was unable to reinstall {uninstallable_components[0]}. Removing 'modules.json' entry")
            else:
                log.info(
                    f"Was unable to reinstall some {component_type}. Removing 'modules.json' entries: {', '.join(uninstallable_components)}"
                )
            if self.modules_json is None:
                raise UserWarning("No modules.json file found")
            for (repo_url, install_dir), component_entries in remove_from_mod_json.items():
                for component in component_entries:
                    self.modules_json["repos"][repo_url][component_type][install_dir].pop(component)
                if len(self.modules_json["repos"][repo_url][component_type][install_dir]) == 0:
                    self.modules_json["repos"].pop(repo_url)

    def resolve_missing_from_modules_json(self, missing_from_modules_json, component_type):
        format_missing = [f"'{dir}'" for dir in missing_from_modules_json]
        if len(format_missing) == 1:
            log.info(
                f"Recomputing commit SHA for {component_type[:-1]} {format_missing[0]} which was missing from 'modules.json'"
            )
        else:
            log.info(
                f"Recomputing commit SHAs for {component_type} which were missing from 'modules.json': {', '.join(format_missing)}"
            )
        assert self.modules_json is not None  # mypy
        # Get the remotes we are missing
        tracked_repos = {repo_url: (repo_entry) for repo_url, repo_entry in self.modules_json["repos"].items()}
        repos, _ = self.get_pipeline_module_repositories(component_type, self.modules_dir, tracked_repos)

        # Get tuples of components that miss installation and their install directory

        def components_with_repos():
            for directory in missing_from_modules_json:
                for repo_url in repos:
                    modules_repo = ModulesRepo(repo_url)
                    paths_in_directory = []
                    repo_url_path = Path(
                        self.modules_dir,
                        modules_repo.repo_path,
                    )
                    for dir_name, _, _ in os.walk(repo_url_path):
                        if component_type == "modules":
                            if len(Path(directory).parts) > 1:  # The module name is TOOL/SUBTOOL
                                paths_in_directory.append(str(Path(*Path(dir_name).parts[-2:])))
                                pass
                        paths_in_directory.append(Path(dir_name).parts[-1])
                    if directory in paths_in_directory:
                        yield (modules_repo.repo_path, directory)

        # Add all components into a dictionary with install directories
        repos_with_components = {}
        for install_dir, component in components_with_repos():
            if install_dir not in repos_with_components:
                repos_with_components[install_dir] = []
            repos_with_components[install_dir].append(component)

        for install_dir, components in repos_with_components.items():
            remote_url = [
                url
                for url, content in repos.items()
                for comp_type, install_directories in content.items()
                if install_dir in install_directories
            ][0]
            repo_entry = self.determine_branches_and_shas(component_type, install_dir, remote_url, components)
            try:
                self.modules_json["repos"][remote_url][component_type][install_dir].update(repo_entry)
            except KeyError:
                try:
                    self.modules_json["repos"][remote_url][component_type].update({install_dir: repo_entry})
                except KeyError:
                    try:
                        self.modules_json["repos"][remote_url].update(
                            {
                                component_type: {
                                    install_dir: repo_entry,
                                }
                            }
                        )
                    except KeyError:
                        self.modules_json["repos"].update(
                            {
                                remote_url: {
                                    component_type: {
                                        install_dir: repo_entry,
                                    }
                                }
                            }
                        )

    def recreate_dependencies(self, repo, org, subworkflow):
        """
        Try to recreate the installed_by entries for subworkflows.
        Remove self installation entry from dependencies, assuming that the modules.json has been freshly created,
        i.e., no module or subworkflow has been installed by the user in the meantime
        """

        sw_name = subworkflow["name"]
        sw_path = Path(self.subworkflows_dir, org, sw_name)
        dep_mods, dep_subwfs = get_components_to_install(sw_path)
        assert self.modules_json is not None  # mypy
        for dep_mod in dep_mods:
            name = dep_mod["name"]
            current_repo = dep_mod.get("git_remote", repo)
            current_org = dep_mod.get("org_path", org)
            installed_by = self.modules_json["repos"][current_repo]["modules"][current_org][name]["installed_by"]
            if installed_by == ["modules"]:
                self.modules_json["repos"][repo]["modules"][org][dep_mod]["installed_by"] = []
            if sw_name not in installed_by:
                self.modules_json["repos"][repo]["modules"][org][dep_mod]["installed_by"].append(sw_name)

        for dep_subwf in dep_subwfs:
            name = dep_subwf["name"]
            current_repo = dep_subwf.get("git_remote", repo)
            current_org = dep_subwf.get("org_path", org)
            installed_by = self.modules_json["repos"][current_repo]["subworkflows"][current_org][name]["installed_by"]
            if installed_by == ["subworkflows"]:
                self.modules_json["repos"][repo]["subworkflows"][org][dep_subwf]["installed_by"] = []
            if sw_name not in installed_by:
                self.modules_json["repos"][repo]["subworkflows"][org][dep_subwf]["installed_by"].append(sw_name)
            self.recreate_dependencies(repo, org, dep_subwf)
