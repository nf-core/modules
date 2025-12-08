import logging
import os
import re
import shutil
from pathlib import Path

import git
import rich
from git.exc import GitCommandError, InvalidGitRepositoryError
from packaging.version import Version

import nf_core
import nf_core.modules.modules_utils
from nf_core.pipelines.download.utils import DownloadError
from nf_core.synced_repo import RemoteProgressbar, SyncedRepo
from nf_core.utils import (
    NFCORE_CACHE_DIR,
    NFCORE_DIR,
)

log = logging.getLogger(__name__)


class WorkflowRepo(SyncedRepo):
    """
    An object to store details about a locally cached workflow repository.

    Important Attributes:
        fullname: The full name of the repository, ``nf-core/{self.pipelinename}``.
        local_repo_dir (str): The local directory, where the workflow is cloned into. Defaults to ``$HOME/.cache/nf-core/nf-core/{self.pipeline}``.

    """

    def __init__(
        self,
        remote_url,
        revision,
        commit,
        additional_tags,
        location=None,
        hide_progress=False,
        in_cache=True,
    ):
        """
        Initializes the object and clones the workflows git repository if it is not already present

        Args:
            remote_url (str): The URL of the remote repository. Defaults to None.
            self.revision (list of str): The revisions to include. A list of strings.
            commits (dict of str): The checksums to linked with the revisions.
            no_pull (bool, optional): Whether to skip the pull step. Defaults to False.
            hide_progress (bool, optional): Whether to hide the progress bar. Defaults to False.
            in_cache (bool, optional): Whether to clone the repository from the cache. Defaults to False.
        """
        self.remote_url = remote_url
        if isinstance(revision, str):
            self.revision = [revision]
        elif isinstance(revision, list):
            self.revision = [*revision]
        else:
            self.revision = []
        if isinstance(commit, str):
            self.commit = [commit]
        elif isinstance(commit, list):
            self.commit = [*commit]
        else:
            self.commit = []
        self.fullname = nf_core.modules.modules_utils.repo_full_name_from_remote(self.remote_url)
        self.retries = 0  # retries for setting up the locally cached repository
        self.hide_progress = hide_progress

        self.setup_local_repo(remote=remote_url, location=location, in_cache=in_cache)

        # additional tags to be added to the repository
        self.additional_tags = additional_tags if additional_tags else None

    def __repr__(self):
        """Called by print, creates representation of object"""
        return f"<Locally cached repository: {self.fullname}, revisions {', '.join(self.revision)}\n cached at: {self.local_repo_dir}>"

    @property
    def heads(self):
        return self.repo.heads

    @property
    def tags(self):
        return self.repo.tags

    def access(self):
        if self.local_repo_dir.exists():
            return self.local_repo_dir
        else:
            return None

    def checkout(self, commit):
        return super().checkout(commit)

    def get_remote_branches(self, remote_url):
        return super().get_remote_branches(remote_url)

    def retry_setup_local_repo(self, skip_confirm=False):
        self.retries += 1
        if skip_confirm or rich.prompt.Confirm.ask(
            f"[violet]Delete local cache '{self.local_repo_dir}' and try again?"
        ):
            if (
                self.retries > 1
            ):  # One unconfirmed retry is acceptable, but prevent infinite loops without user interaction.
                raise DownloadError(
                    f"Errors with locally cached repository of '{self.fullname}'. Please delete '{self.local_repo_dir}' manually and try again."
                )
            if not skip_confirm:  # Feedback to user for manual confirmation.
                log.info(f"Removing '{self.local_repo_dir}'")
            shutil.rmtree(self.local_repo_dir)
            self.setup_local_repo(self.remote_url, in_cache=False)
        else:
            raise DownloadError("Exiting due to error with locally cached Git repository.")

    def setup_local_repo(self, remote, location=None, in_cache=True):
        """
        Sets up the local git repository. If the repository has been cloned previously, it
        returns a git.Repo object of that clone. Otherwise it tries to clone the repository from
        the provided remote URL and returns a git.Repo of the new clone.

        Args:
            remote (str): git url of remote
            location (Path): location where the clone should be created/cached.
            in_cache (bool, optional): Whether to clone the repository from the cache. Defaults to False.
        Sets self.repo
        """
        if location:
            self.local_repo_dir = location / self.fullname
        else:
            self.local_repo_dir = Path(NFCORE_DIR) if not in_cache else Path(NFCORE_CACHE_DIR, self.fullname)

        try:
            if not self.local_repo_dir.exists():
                try:
                    pbar = rich.progress.Progress(
                        "[bold blue]{task.description}",
                        rich.progress.BarColumn(bar_width=None),
                        "[bold yellow]{task.fields[state]}",
                        transient=True,
                        disable=os.environ.get("HIDE_PROGRESS", None) is not None or self.hide_progress,
                    )
                    with pbar:
                        self.repo = git.Repo.clone_from(
                            remote,
                            self.local_repo_dir,
                            progress=RemoteProgressbar(pbar, self.fullname, self.remote_url, "Cloning"),
                        )
                    super().update_local_repo_status(self.fullname, True)
                except GitCommandError:
                    raise DownloadError(f"Failed to clone from the remote: `{remote}`")
            else:
                self.repo = git.Repo(self.local_repo_dir)

                if super().no_pull_global:
                    super().update_local_repo_status(self.fullname, True)
                # If the repo is already cloned, fetch the latest changes from the remote
                if not super().local_repo_synced(self.fullname):
                    pbar = rich.progress.Progress(
                        "[bold blue]{task.description}",
                        rich.progress.BarColumn(bar_width=None),
                        "[bold yellow]{task.fields[state]}",
                        transient=True,
                        disable=os.environ.get("HIDE_PROGRESS", None) is not None or self.hide_progress,
                    )
                    with pbar:
                        self.repo.remotes.origin.fetch(
                            progress=RemoteProgressbar(pbar, self.fullname, self.remote_url, "Pulling")
                        )
                    super().update_local_repo_status(self.fullname, True)

        except (GitCommandError, InvalidGitRepositoryError) as e:
            log.error(f"[red]Could not set up local cache of modules repository:[/]\n{e}\n")
            self.retry_setup_local_repo()

    def tidy_tags_and_branches(self):
        """
        Function to delete all tags and branches that are not of interest to the downloader.
        This allows a clutter-free experience in Seqera Platform. The untagged commits are evidently still available.

        However, due to local caching, the downloader might also want access to revisions that had been deleted before.
        In that case, don't bother with re-adding the tags and rather download  anew from Github.
        """
        if self.revision and self.repo and self.repo.tags:
            # create a set to keep track of the revisions to process & check
            desired_revisions = set(self.revision)

            # determine what needs pruning
            tags_to_remove = {tag for tag in self.repo.tags if tag.name not in desired_revisions.union({"latest"})}
            heads_to_remove = {head for head in self.repo.heads if head.name not in desired_revisions.union({"latest"})}

            try:
                # delete unwanted tags from repository
                for tag in tags_to_remove:
                    self.repo.delete_tag(tag)

                # switch to a revision that should be kept, because deleting heads fails, if they are checked out (e.g. "main")
                self.checkout(self.revision[0])

                # delete unwanted heads/branches from repository
                for head in heads_to_remove:
                    self.repo.delete_head(head)

                # ensure all desired revisions/branches are available
                for revision in desired_revisions:
                    if not self.repo.is_valid_object(revision):
                        self.checkout(revision)
                        self.repo.create_head(revision, revision)
                        if self.repo.head.is_detached:
                            self.repo.head.reset(index=True, working_tree=True)

                # no branch exists, but one is required for Seqera Platform's UI to display revisions correctly). Thus, "latest" will be created.
                if not bool(self.repo.heads):
                    if self.repo.is_valid_object("latest"):
                        # "latest" exists as tag but not as branch
                        self.repo.create_head("latest", "latest")  # create a new head for latest
                        self.checkout("latest")
                    else:
                        # desired revisions may contain arbitrary branch names that do not correspond to valid semantic versioning patterns.
                        valid_versions = [
                            Version(v) for v in desired_revisions if re.match(r"\d+\.\d+(?:\.\d+)*(?:[\w\-_])*", v)
                        ]
                        # valid versions sorted in ascending order, last will be aliased as "latest".
                        latest = sorted(valid_versions)[-1]
                        self.repo.create_head("latest", str(latest))
                        self.checkout(latest)
                    if self.repo.head.is_detached:
                        self.repo.head.reset(index=True, working_tree=True)

                # Apply the custom additional tags to the repository
                self.__add_additional_tags()

                # get all tags and available remote_branches
                completed_revisions = {revision.name for revision in self.repo.heads + self.repo.tags}

                # verify that all requested revisions are available.
                # a local cache might lack revisions that were deleted during a less comprehensive previous download.
                if bool(desired_revisions - completed_revisions):
                    log.info(
                        f"Locally cached version of the pipeline lacks selected revisions {', '.join(desired_revisions - completed_revisions)}. Downloading anew from GitHub..."
                    )
                    self.retry_setup_local_repo(skip_confirm=True)
                    self.tidy_tags_and_branches()
            except (GitCommandError, InvalidGitRepositoryError) as e:
                log.error(f"[red]Adapting your pipeline download unfortunately failed:[/]\n{e}\n")
                self.retry_setup_local_repo(skip_confirm=True)
                raise DownloadError(e) from e

    # "Private" method to add the additional custom tags to the repository.
    def __add_additional_tags(self) -> None:
        if self.additional_tags:
            # example.com is reserved by the Internet Assigned Numbers Authority (IANA)  as special-use domain names for documentation purposes.
            # Although "dev-null" is a syntactically-valid local-part that is equally valid for delivery,
            # and only the receiving MTA can decide whether to accept it, it is to my best knowledge configured with
            # a Postfix discard mail delivery agent (https://www.postfix.org/discard.8.html), so incoming mails should be sinkholed.
            self.ensure_git_user_config(
                f"nf-core pipelines download v{nf_core.__version__}",
                "dev-null@example.com",
            )

            for additional_tag in self.additional_tags:
                # A valid git branch or tag name can contain alphanumeric characters, underscores, hyphens, and dots.
                # But it must not start with a dot, hyphen or underscore and also cannot contain two consecutive dots.
                if re.match(r"^\w[\w_.-]+={1}\w[\w_.-]+$", additional_tag) and ".." not in additional_tag:
                    anchor, tag = additional_tag.split("=")
                    if self.repo.is_valid_object(anchor) and not self.repo.is_valid_object(tag):
                        try:
                            self.repo.create_tag(
                                tag,
                                ref=anchor,
                                message=f"Synonynmous tag to {anchor}; added by `nf-core pipelines download`.",
                            )
                        except (GitCommandError, InvalidGitRepositoryError) as e:
                            log.error(f"[red]Additional tag(s) could not be applied:[/]\n{e}\n")
                    else:
                        if not self.repo.is_valid_object(anchor):
                            log.error(
                                f"[red]Adding tag '{tag}' to '{anchor}' failed.[/]\n Mind that '{anchor}' must be a valid git reference that resolves to a commit."
                            )
                        if self.repo.is_valid_object(tag):
                            log.error(
                                f"[red]Adding tag '{tag}' to '{anchor}' failed.[/]\n Mind that '{tag}' must not exist hitherto."
                            )
                else:
                    log.error(f"[red]Could not apply invalid `--tag` specification[/]: '{additional_tag}'")

    def bare_clone(self, destination: Path):
        if self.repo:
            try:
                destfolder = destination.parent.absolute()
                if not destfolder.exists():
                    destfolder.mkdir()
                if destination.exists():
                    shutil.rmtree(destination)
                self.repo.clone(str(destfolder), bare=True)
            except (OSError, GitCommandError, InvalidGitRepositoryError) as e:
                log.error(f"[red]Failure to create the pipeline download[/]\n{e}\n")
