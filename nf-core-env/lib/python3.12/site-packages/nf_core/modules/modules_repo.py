import logging
import os
import shutil
from pathlib import Path

import git
import rich
import rich.progress
import rich.prompt
from git.exc import GitCommandError, InvalidGitRepositoryError

import nf_core.modules.modules_utils
from nf_core.components.constants import NF_CORE_MODULES_NAME, NF_CORE_MODULES_REMOTE
from nf_core.synced_repo import RemoteProgressbar, SyncedRepo
from nf_core.utils import NFCORE_CACHE_DIR, NFCORE_DIR, load_tools_config

log = logging.getLogger(__name__)


class ModulesRepo(SyncedRepo):
    """
    An object to store details about the repository being used for modules.

    Used by the `nf-core modules` top-level command with -r and -b flags,
    so that this can be used in the same way by all sub-commands.

    We keep track of the pull-status of the different installed repos in
    the static variable local_repo_status. This is so we don't need to
    pull a remote several times in one command.
    """

    local_repo_statuses = {}
    no_pull_global = False

    def __init__(
        self,
        remote_url: str | None = None,
        branch: str | None = None,
        no_pull: bool = False,
        hide_progress: bool = False,
    ) -> None:
        """
        Initializes the object and clones the git repository if it is not already present
        """

        # This allows us to set this one time and then keep track of the user's choice
        ModulesRepo.no_pull_global |= no_pull

        # Check if the remote seems to be well formed
        if remote_url is None:
            remote_url = NF_CORE_MODULES_REMOTE

        self.remote_url = remote_url

        self.fullname = nf_core.modules.modules_utils.repo_full_name_from_remote(self.remote_url)

        self.setup_local_repo(remote_url, branch, hide_progress)

        config_fn, repo_config = load_tools_config(self.local_repo_dir)
        if config_fn is None or repo_config is None:
            raise UserWarning(f"Could not find a configuration file in {self.local_repo_dir}")
        try:
            self.repo_path = repo_config.org_path
        except KeyError:
            raise UserWarning(f"'org_path' key not present in {config_fn.name}")

        # Verify that the repo seems to be correctly configured
        if self.repo_path != NF_CORE_MODULES_NAME or self.branch:
            self.verify_branch()
        if self.repo_path is None:
            raise UserWarning(f"Could not find the org_path in the configuration file: {config_fn.name}")
        # Convenience variable
        self.modules_dir = Path(self.local_repo_dir, "modules", self.repo_path)
        self.subworkflows_dir = Path(self.local_repo_dir, "subworkflows", self.repo_path)

        self.avail_module_names = None

    def gitless_repo(self):
        gitless_repo_url = self.remote_url
        if self.remote_url and ".git" in self.remote_url:
            gitless_repo_url = gitless_repo_url[:-4]
        return gitless_repo_url

    def setup_local_repo(self, remote, branch, hide_progress=True, in_cache=False):
        """
        Sets up the local git repository. If the repository has been cloned previously, it
        returns a git.Repo object of that clone. Otherwise it tries to clone the repository from
        the provided remote URL and returns a git.Repo of the new clone.

        Args:
            remote (str): git url of remote
            branch (str): name of branch to use
        Sets self.repo
        """
        self.local_repo_dir = Path(NFCORE_DIR if not in_cache else NFCORE_CACHE_DIR, self.fullname)
        try:
            if not os.path.exists(self.local_repo_dir):
                try:
                    pbar = rich.progress.Progress(
                        "[bold blue]{task.description}",
                        rich.progress.BarColumn(bar_width=None),
                        "[bold yellow]{task.fields[state]}",
                        transient=True,
                        disable=hide_progress or os.environ.get("HIDE_PROGRESS", None) is not None,
                    )
                    with pbar:
                        self.repo = git.Repo.clone_from(
                            remote,
                            self.local_repo_dir,
                            progress=RemoteProgressbar(pbar, self.fullname, self.remote_url, "Cloning"),
                        )
                    ModulesRepo.update_local_repo_status(self.fullname, True)
                except GitCommandError:
                    raise LookupError(f"Failed to clone from the remote: `{remote}`")
                # Verify that the requested branch exists by checking it out
                self.setup_branch(branch)
            else:
                self.repo = git.Repo(self.local_repo_dir)

                if ModulesRepo.no_pull_global:
                    ModulesRepo.update_local_repo_status(self.fullname, True)
                # If the repo is already cloned, fetch the latest changes from the remote
                if not ModulesRepo.local_repo_synced(self.fullname):
                    pbar = rich.progress.Progress(
                        "[bold blue]{task.description}",
                        rich.progress.BarColumn(bar_width=None),
                        "[bold yellow]{task.fields[state]}",
                        transient=True,
                        disable=hide_progress or os.environ.get("HIDE_PROGRESS", None) is not None,
                    )
                    with pbar:
                        self.repo.remotes.origin.fetch(
                            progress=RemoteProgressbar(pbar, self.fullname, self.remote_url, "Pulling")
                        )
                    ModulesRepo.update_local_repo_status(self.fullname, True)

                # Before verifying the branch, fetch the changes
                # Verify that the requested branch exists by checking it out
                self.setup_branch(branch)

                # Now merge the changes
                tracking_branch = self.repo.active_branch.tracking_branch()
                if tracking_branch is None:
                    raise LookupError(f"There is no remote tracking branch '{self.branch}' in '{self.remote_url}'")
                self.repo.git.merge(tracking_branch.name)
        except (GitCommandError, InvalidGitRepositoryError) as e:
            log.error(f"[red]Could not set up local cache of modules repository:[/]\n{e}\n")
            if rich.prompt.Confirm.ask(f"[violet]Delete local cache '{self.local_repo_dir}' and try again?"):
                log.info(f"Removing '{self.local_repo_dir}'")
                shutil.rmtree(self.local_repo_dir)
                self.setup_local_repo(remote, branch, hide_progress)
            else:
                raise LookupError("Exiting due to error with local modules git repo")
