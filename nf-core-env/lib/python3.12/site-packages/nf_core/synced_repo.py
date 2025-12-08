import filecmp
import logging
import os
import shutil
from collections.abc import Iterable
from configparser import NoOptionError, NoSectionError
from pathlib import Path

import git
from git.exc import GitCommandError

from nf_core.components.constants import NF_CORE_MODULES_DEFAULT_BRANCH, NF_CORE_MODULES_NAME, NF_CORE_MODULES_REMOTE
from nf_core.utils import load_tools_config

log = logging.getLogger(__name__)


class RemoteProgressbar(git.RemoteProgress):
    """
    An object to create a progressbar for when doing an operation with the remote.
    Note that an initialized rich Progress (progress bar) object must be passed
    during initialization.
    """

    def __init__(self, progress_bar, repo_name, remote_url, operation):
        """
        Initializes the object and adds a task to the progressbar passed as 'progress_bar'

        Args:
            progress_bar (rich.progress.Progress): A rich progress bar object
            repo_name (str): Name of the repository the operation is performed on
            remote_url (str): Git URL of the repository the operation is performed on
            operation (str): The operation performed on the repository, i.e. 'Pulling', 'Cloning' etc.
        """
        super().__init__()
        self.progress_bar = progress_bar
        self.tid = self.progress_bar.add_task(
            f"{operation} from [bold green]'{repo_name}'[/bold green] ([link={remote_url}]{remote_url}[/link])",
            start=False,
            state="Waiting for response",
        )

    def update(self, op_code, cur_count, max_count=None, message=""):
        """
        Overrides git.RemoteProgress.update.
        Called every time there is a change in the remote operation
        """
        if not self.progress_bar.tasks[self.tid].started:
            self.progress_bar.start_task(self.tid)
        if cur_count is not None and max_count is not None:
            cur_count = float(cur_count)
            max_count = float(max_count)
            state = f"{cur_count / max_count * 100:.1f}%"
        else:
            state = "Unknown"

        self.progress_bar.update(self.tid, total=max_count, completed=cur_count, state=state)


class SyncedRepo:
    """
    An object to store details about a locally cached code repository.
    """

    local_repo_statuses: dict[str, bool] = {}
    no_pull_global = False

    @staticmethod
    def local_repo_synced(repo_name):
        """
        Checks whether a local repo has been cloned/pull in the current session
        """
        return SyncedRepo.local_repo_statuses.get(repo_name, False)

    @staticmethod
    def update_local_repo_status(repo_name, up_to_date):
        """
        Updates the clone/pull status of a local repo
        """
        SyncedRepo.local_repo_statuses[repo_name] = up_to_date

    @staticmethod
    def get_remote_branches(remote_url):
        """
        Get all branches from a remote repository

        Args:
            remote_url (str): The git url to the remote repository

        Returns:
            (set[str]): All branches found in the remote
        """
        try:
            unparsed_branches = git.Git().ls_remote(remote_url)
        except git.GitCommandError:
            raise LookupError(f"Was unable to fetch branches from '{remote_url}'")
        else:
            branches = {}
            for branch_info in unparsed_branches.split("\n"):
                sha, name = branch_info.split("\t")
                if name != "HEAD":
                    # The remote branches are shown as 'ref/head/branch'
                    branch_name = Path(name).stem
                    branches[sha] = branch_name
            return set(branches.values())

    def __init__(self, remote_url=None, branch=None, no_pull=False, hide_progress=False):
        """
        Initializes the object and clones the git repository if it is not already present
        """

        # This allows us to set this one time and then keep track of the user's choice
        SyncedRepo.no_pull_global |= no_pull

        # Check if the remote seems to be well formed
        if remote_url is None:
            remote_url = NF_CORE_MODULES_REMOTE

        self.remote_url = remote_url
        self.fullname = None
        self.local_repo_dir = None

        self.repo = None
        # TODO: SyncedRepo doesn't have this method and both the ModulesRepo and
        # the WorkflowRepo define their own including custom init methods. This needs
        # fixing.
        self.setup_local_repo(remote_url, branch, hide_progress)

        if self.local_repo_dir is None:
            raise ValueError("Repository not initialized")
        else:
            config_fn, repo_config = load_tools_config(self.local_repo_dir)
            if config_fn is not None and repo_config is not None:
                try:
                    self.repo_path = repo_config.org_path
                except KeyError:
                    raise UserWarning(f"'org_path' key not present in {config_fn.name}")

            # Verify that the repo seems to be correctly configured
            if self.repo_path != NF_CORE_MODULES_NAME or self.branch:
                self.verify_branch()

            # Convenience variable
            self.modules_dir = Path(self.local_repo_dir, "modules", self.repo_path)
            self.subworkflows_dir = Path(self.local_repo_dir, "subworkflows", self.repo_path)

            self.avail_module_names = None

    def __repr__(self) -> str:
        return f"SyncedRepo({self.remote_url}, {self.branch})"

    def setup_local_repo(self, remote_url, branch, hide_progress):
        pass

    def verify_sha(self, prompt, sha):
        """
        Verify that 'sha' and 'prompt' arguments are not provided together.
        Verify that the provided SHA exists in the repo.

        Arguments:
            prompt (bool):              prompt asking for SHA
            sha (str):                  provided sha
        """
        if prompt and sha is not None:
            log.error("Cannot use '--sha' and '--prompt' at the same time!")
            return False

        if sha:
            if not self.sha_exists_on_branch(sha):
                log.error(f"Commit SHA '{sha}' doesn't exist in '{self.remote_url}'")
                return False

        return True

    def setup_branch(self, branch):
        """
        Verify that we have a branch and otherwise use the default one.
        The branch is then checked out to verify that it exists in the repo.

        Args:
            branch (str): Name of branch
        """
        if branch is None:
            # Don't bother fetching default branch if we're using nf-core
            if self.remote_url == NF_CORE_MODULES_REMOTE:
                self.branch = NF_CORE_MODULES_DEFAULT_BRANCH
            else:
                self.branch = self.get_default_branch()
        else:
            self.branch = branch

        # Verify that the branch exists by checking it out
        self.branch_exists()

    def get_default_branch(self):
        """
        Gets the default branch for the repo (the branch origin/HEAD is pointing to)
        """
        origin_head = next(ref for ref in self.repo.refs if ref.name == "origin/HEAD")
        _, branch = origin_head.ref.name.split("/")
        return branch

    def branch_exists(self):
        """
        Verifies that the branch exists in the repository by trying to check it out
        """
        try:
            self.checkout_branch()
        except GitCommandError:
            raise LookupError(f"Branch '{self.branch}' not found in '{self.remote_url}'")

    def verify_branch(self):
        """
        Verifies the active branch conforms to the correct directory structure
        """
        dir_names = os.listdir(self.local_repo_dir)
        if "modules" not in dir_names:
            err_str = f"Repository '{self.remote_url}' ({self.branch}) does not contain the 'modules/' directory"
            if "software" in dir_names:
                err_str += (
                    ".\nAs of nf-core/tools version 2.0, the 'software/' directory should be renamed to 'modules/'"
                )
            raise LookupError(err_str)

    def checkout_branch(self):
        """
        Checks out the specified branch of the repository
        """
        try:
            self.repo.git.checkout(self.branch)
        except GitCommandError as e:
            if (
                self.fullname
                and "modules" in self.fullname
                and "Your local changes to the following files would be overwritten by checkout" in str(e)
            ):
                log.debug(f"Overwriting local changes in '{self.local_repo_dir}'")
                self.repo.git.checkout(self.branch, force=True)
            else:
                raise e

    def checkout(self, commit):
        """
        Checks out the repository at the requested commit

        Args:
            commit (str): Git SHA of the commit
        """
        try:
            self.repo.git.checkout(commit)
        except GitCommandError as e:
            if (
                self.fullname
                and "modules" in self.fullname
                and "Your local changes to the following files would be overwritten by checkout" in str(e)
            ):
                log.debug(f"Overwriting local changes in '{self.local_repo_dir}'")
                self.repo.git.checkout(self.branch, force=True)
            else:
                raise e

    def component_exists(self, component_name, component_type, checkout=True, commit=None):
        """
        Check if a module/subworkflow exists in the branch of the repo

        Args:
            component_name (str): The name of the module/subworkflow

        Returns:
            (bool): Whether the module/subworkflow exists in this branch of the repository
        """
        return component_name in self.get_avail_components(component_type, checkout=checkout, commit=commit)

    def get_component_dir(self, component_name: str, component_type: str) -> Path:
        """
        Returns the file path of a module/subworkflow directory in the repo.
        Does not verify that the path exists.
        Args:
            component_name (str): The name of the module/subworkflow

        Returns:
            component_path (str): The path of the module/subworkflow in the local copy of the repository
        """
        if component_type == "modules":
            return Path(self.modules_dir, component_name)
        elif component_type == "subworkflows":
            return Path(self.subworkflows_dir, component_name)
        else:
            raise ValueError(f"Invalid component type: {component_type}")

    def install_component(self, component_name: str, install_dir: str | Path, commit: str, component_type: str) -> bool:
        """
        Install the module/subworkflow files into a pipeline at the given commit

        Args:
            component_name (str): The name of the module/subworkflow
            install_dir (str): The path where the module/subworkflow should be installed
            commit (str): The git SHA for the version of the module/subworkflow to be installed
            component_type (str): Either 'modules' or 'subworkflows'

        Returns:
            (bool): Whether the operation was successful or not
        """
        # Check out the repository at the requested ref
        try:
            self.checkout(commit)
        except git.GitCommandError:
            return False

        # Check if the module/subworkflow exists in the branch
        if not self.component_exists(component_name, component_type, checkout=False):
            log.error(
                f"The requested {component_type[:-1]} does not exists in the branch '{self.branch}' of {self.remote_url}'"
            )
            return False

        # Copy the files from the repo to the install folder
        shutil.copytree(self.get_component_dir(component_name, component_type), Path(install_dir, component_name))

        # Switch back to the tip of the branch
        self.checkout_branch()
        return True

    def component_files_identical(self, component_name, base_path, commit, component_type):
        """
        Checks whether the module or subworkflow files in a pipeline are identical to the ones in the remote
        Args:
            component_name (str): The name of the module or subworkflow
            base_path (str): The path to the module/subworkflow in the pipeline

        Returns:
            (bool): Whether the pipeline files are identical to the repo files
        """
        if commit is None:
            self.checkout_branch()
        else:
            self.checkout(commit)
        component_files = ["main.nf", "meta.yml"]
        files_identical = {file: True for file in component_files}
        component_dir = self.get_component_dir(component_name, component_type)
        for file in component_files:
            try:
                files_identical[file] = filecmp.cmp(Path(component_dir, file), Path(base_path, file))
            except FileNotFoundError:
                log.debug(f"Could not open file: {Path(component_dir, file)}")
                continue
        self.checkout_branch()
        return files_identical

    def ensure_git_user_config(self, default_name: str, default_email: str) -> None:
        if self.repo is None:
            raise ValueError("Repository not initialized")
        try:
            with self.repo.config_reader() as git_config:
                user_name = git_config.get_value("user", "name", default=None)
                user_email = git_config.get_value("user", "email", default=None)
        except (NoOptionError, NoSectionError):
            user_name = user_email = None

        if not user_name or not user_email:
            with self.repo.config_writer() as git_config:
                if not user_name:
                    git_config.set_value("user", "name", default_name)
                if not user_email:
                    git_config.set_value("user", "email", default_email)

    def get_component_git_log(
        self, component_name: str | Path, component_type: str, depth: int | None = None
    ) -> Iterable[dict[str, str]]:
        """
        Fetches the commit history the of requested module/subworkflow since a given date. The default value is
        not arbitrary - it is the last time the structure of the nf-core/modules repository was had an
        update breaking backwards compatibility.
        Args:
            component_name (str): Name of module/subworkflow
            modules_repo (SyncedRepo): A SyncedRepo object configured for the repository in question

        Returns:
            ( dict ): Iterator of commit SHAs and associated (truncated) message
        """
        if self.repo is None:
            raise ValueError("Repository not initialized")
        self.checkout_branch()
        component_path = Path(component_type, self.repo_path, component_name)

        commits_new_iter = self.repo.iter_commits(max_count=depth, paths=component_path)
        commits_old_iter = []
        if component_type == "modules":
            # Grab commits also from previous modules structure
            old_component_path = Path("modules", component_name)
            commits_old_iter = self.repo.iter_commits(max_count=depth, paths=old_component_path)

        try:
            commits_old = [{"git_sha": commit.hexsha, "trunc_message": commit.message} for commit in commits_old_iter]
            commits_new = [{"git_sha": commit.hexsha, "trunc_message": commit.message} for commit in commits_new_iter]
        except git.GitCommandError as e:
            log.error(
                f"Git error: {e}\n"
                "To solve this, you can try to remove the cloned rempository and run the command again.\n"
                f"This repository is typically found at `{self.local_repo_dir}`"
            )
            raise UserWarning
        commits = iter(commits_new + commits_old)

        return commits

    def get_latest_component_version(self, component_name, component_type):
        """
        Returns the latest commit in the repository
        """
        try:
            git_logs = list(self.get_component_git_log(component_name, component_type, depth=1))
            if not git_logs:
                return None
            return git_logs[0]["git_sha"]
        except Exception as e:
            log.debug(f"Could not get latest version of {component_name}: {e}")
            return None

    def sha_exists_on_branch(self, sha):
        """
        Verifies that a given commit sha exists on the branch
        """
        if self.repo is None:
            raise ValueError("Repository not initialized")
        self.checkout_branch()
        return sha in (commit.hexsha for commit in self.repo.iter_commits())

    def get_commit_info(self, sha):
        """
        Fetches metadata about the commit (dates, message, etc.)
        Args:
            commit_sha (str): The SHA of the requested commit
        Returns:
            message (str): The commit message for the requested commit
            date (str): The commit date for the requested commit
        Raises:
            LookupError: If the search for the commit fails
        """
        if self.repo is None:
            raise ValueError("Repository not initialized")
        self.checkout_branch()
        for commit in self.repo.iter_commits():
            if commit.hexsha == sha:
                message = commit.message.splitlines()[0]
                date_obj = commit.committed_datetime
                date = str(date_obj.date())
                return message, date
        raise LookupError(f"Commit '{sha}' not found in the '{self.remote_url}'")

    def get_avail_components(self, component_type: str, checkout: bool = True, commit: str | None = None) -> list[str]:
        """
        Gets the names of the modules/subworkflows in the repository. They are detected by
        checking which directories have a 'main.nf' file

        Returns:
            ([ str ]): The module/subworkflow names
        """
        if checkout:
            self.checkout_branch()
        if commit is not None:
            self.checkout(commit)
        # Get directory
        if component_type == "modules":
            directory = self.modules_dir
        elif component_type == "subworkflows":
            directory = self.subworkflows_dir
        # Module/Subworkflow directories are characterized by having a 'main.nf' file
        avail_component_names = [
            str(Path(dirpath).relative_to(directory)) for dirpath, _, files in os.walk(directory) if "main.nf" in files
        ]
        return avail_component_names

    def get_meta_yml(self, component_type, module_name):
        """
        Returns the contents of the 'meta.yml' file of a module

        Args:
            module_name (str): The name of the module

        Returns:
            (str): The contents of the file in text format
        """
        self.checkout_branch()
        if component_type == "modules":
            path = Path(self.modules_dir, module_name, "meta.yml")
        elif component_type == "subworkflows":
            path = Path(self.subworkflows_dir, module_name, "meta.yml")
        else:
            raise ValueError(f"Invalid component type: {component_type}")
        if not path.exists():
            return None
        with open(path) as fh:
            contents = fh.read()
        return contents
