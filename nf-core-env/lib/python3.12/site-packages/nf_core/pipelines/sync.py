"""Synchronise a pipeline TEMPLATE branch with the template."""

import json
import logging
import os
import re
from pathlib import Path
from typing import Any

import git
import questionary
import requests
import requests.auth
import requests_cache
import rich
import yaml
from git import GitCommandError, InvalidGitRepositoryError

import nf_core
import nf_core.pipelines.create.create
import nf_core.pipelines.list
import nf_core.utils
from nf_core.pipelines.lint_utils import dump_yaml_with_prettier

log = logging.getLogger(__name__)


class SyncExceptionError(Exception):
    """Exception raised when there was an error with TEMPLATE branch synchronisation"""

    pass


class PullRequestExceptionError(Exception):
    """Exception raised when there was an error creating a Pull-Request on GitHub.com"""

    pass


class PipelineSync:
    """Object to hold syncing information and results.

    Args:
        pipeline_dir (str): The path to the Nextflow pipeline root directory
        from_branch (str): The branch to use to fetch config vars. If not set, will use current active branch
        make_pr (bool): Set this to `True` to create a GitHub pull-request with the changes
        gh_username (str): GitHub username
        gh_repo (str): GitHub repository name
        template_yaml_path (str): Path to template.yml file for pipeline creation settings. DEPRECATED
        force_pr (bool): Force the creation of a pull request, even if there are no changes to the template

    Attributes:
        pipeline_dir (str): Path to target pipeline directory
        from_branch (str): Repo branch to use when collecting workflow variables. Default: active branch.
        original_branch (str): Repo branch that was checked out before we started.
        made_changes (bool): Whether making the new template pipeline introduced any changes
        make_pr (bool): Whether to try to automatically make a PR on GitHub.com
        required_config_vars (list): List of nextflow variables required to make template pipeline
        gh_username (str): GitHub username
        gh_repo (str): GitHub repository name
    """

    def __init__(
        self,
        pipeline_dir: str | Path,
        from_branch: str | None = None,
        make_pr: bool = False,
        gh_repo: str | None = None,
        gh_username: str | None = None,
        template_yaml_path: str | None = None,
        force_pr: bool = False,
        blog_post: str = "",
    ):
        """Initialise syncing object"""

        self.pipeline_dir: Path = Path(pipeline_dir).resolve()
        self.from_branch = from_branch
        self.original_branch = None
        self.original_merge_branch = f"nf-core-template-merge-{nf_core.__version__}"
        self.merge_branch = self.original_merge_branch
        self.made_changes = False
        self.make_pr = make_pr
        self.gh_pr_returned_data: dict = {}
        self.required_config_vars = [
            "manifest.name",
            "manifest.description",
            "manifest.version",
            "manifest.contributors",
        ]
        self.force_pr = force_pr

        self.gh_username = gh_username
        self.gh_repo = gh_repo
        self.pr_url = ""
        self.blog_post = blog_post

        self.config_yml_path, self.config_yml = nf_core.utils.load_tools_config(self.pipeline_dir)
        assert self.config_yml_path is not None and self.config_yml is not None  # mypy
        # Throw deprecation warning if template_yaml_path is set
        if template_yaml_path is not None:
            log.warning(
                f"The `template_yaml_path` argument is deprecated. Saving pipeline creation settings in .nf-core.yml instead. Please remove {template_yaml_path} file."
            )
            if getattr(self.config_yml, "template", None) is not None:
                overwrite_template = questionary.confirm(
                    f"A template section already exists in '{self.config_yml_path}'. Do you want to overwrite?",
                    style=nf_core.utils.nfcore_question_style,
                    default=False,
                ).unsafe_ask()
            if overwrite_template or getattr(self.config_yml, "template", None) is None:
                with open(template_yaml_path) as f:
                    self.config_yml.template = yaml.safe_load(f)
                with open(self.config_yml_path, "w") as fh:
                    yaml.safe_dump(self.config_yml.model_dump(exclude_none=True), fh)
                log.info(f"Saved pipeline creation settings to '{self.config_yml_path}'")
                raise SystemExit(
                    f"Please commit your changes and delete the {template_yaml_path} file. Then run the sync command again."
                )

        # Set up the API auth if supplied on the command line
        self.gh_api = nf_core.utils.gh_api
        self.gh_api.lazy_init()
        if self.gh_username and "GITHUB_AUTH_TOKEN" in os.environ:
            log.debug(f"Authenticating sync as {self.gh_username}")
            self.gh_api.setup_github_auth(
                requests.auth.HTTPBasicAuth(self.gh_username, os.environ["GITHUB_AUTH_TOKEN"])
            )

    def sync(self) -> None:
        """Find workflow attributes, create a new template pipeline on TEMPLATE"""

        # Clear requests_cache so that we don't get stale API responses
        requests_cache.clear()

        log.info(f"Pipeline directory: {self.pipeline_dir}")
        if self.from_branch:
            log.info(f"Using branch '{self.from_branch}' to fetch workflow variables")
        if self.make_pr:
            log.info("Will attempt to automatically create a pull request")

        self.inspect_sync_dir()
        self.get_wf_config()
        self.checkout_template_branch()
        self.delete_tracked_template_branch_files()
        self.make_template_pipeline()
        self.commit_template_changes()

        if not self.made_changes and self.force_pr:
            log.info("No changes made to TEMPLATE, but PR forced")
            self.made_changes = True

        # Push and make a pull request if we've been asked to
        if self.made_changes and self.make_pr or self.force_pr:
            try:
                # Check that we have an API auth token
                if os.environ.get("GITHUB_AUTH_TOKEN", "") == "":
                    raise PullRequestExceptionError("GITHUB_AUTH_TOKEN not set!")

                # Check that we know the github username and repo name
                if self.gh_username is None and self.gh_repo is None:
                    raise PullRequestExceptionError("Could not find GitHub username and repo name")

                self.push_template_branch()
                self.create_merge_base_branch()
                self.push_merge_branch()
                self.make_pull_request()
            except PullRequestExceptionError as e:
                self.reset_target_dir()
                raise PullRequestExceptionError(e)

        self.reset_target_dir()

        if not self.made_changes:
            log.info("No changes made to TEMPLATE - sync complete")
        elif not self.make_pr:
            log.info(
                f"Now try to merge the updates in to your pipeline:\n  cd {self.pipeline_dir}\n  git merge TEMPLATE"
            )

    def inspect_sync_dir(self):
        """Takes a look at the target directory for syncing. Checks that it's a git repo
        and makes sure that there are no uncommitted changes.
        """
        # Check that the pipeline_dir is a git repo
        try:
            self.repo = git.Repo(self.pipeline_dir)
        except InvalidGitRepositoryError:
            raise SyncExceptionError(f"'{self.pipeline_dir}' does not appear to be a git repository")

        # get current branch so we can switch back later
        self.original_branch = self.repo.active_branch.name
        log.info(f"Original pipeline repository branch is '{self.original_branch}'")

        # Check to see if there are uncommitted changes on current branch
        if self.repo.is_dirty(untracked_files=True):
            raise SyncExceptionError(
                "Uncommitted changes found in pipeline directory!\n"
                "Please commit these before running nf-core pipelines sync.\n"
                "(Hint: .gitignored files are ignored.)"
            )

        # Track ignored files to avoid processing them
        self.ignored_files = self._get_ignored_files()

    def get_wf_config(self):
        """Check out the target branch if requested and fetch the nextflow config.
        Check that we have the required config variables.
        """
        # Try to check out target branch (eg. `origin/dev`)
        try:
            if self.from_branch and self.repo.active_branch.name != self.from_branch:
                log.info(f"Checking out workflow branch '{self.from_branch}'")
                self.repo.git.checkout(self.from_branch)
        except GitCommandError:
            raise SyncExceptionError(f"Branch `{self.from_branch}` not found!")

        # If not specified, get the name of the active branch
        if not self.from_branch:
            try:
                self.from_branch = self.repo.active_branch.name
            except GitCommandError as e:
                log.error(f"Could not find active repo branch: {e}")

        # Fetch workflow variables
        log.debug("Fetching workflow config variables")
        self.wf_config = nf_core.utils.fetch_wf_config(Path(self.pipeline_dir))

        # Check that we have the required variables
        for rvar in self.required_config_vars:
            if rvar not in self.wf_config:
                raise SyncExceptionError(f"Workflow config variable `{rvar}` not found!")

    def checkout_template_branch(self):
        """
        Try to check out the origin/TEMPLATE in a new TEMPLATE branch.
        If this fails, try to check out an existing local TEMPLATE branch.
        """
        # Try to check out the `TEMPLATE` branch
        try:
            self.repo.git.checkout("origin/TEMPLATE", b="TEMPLATE")
        except GitCommandError:
            # Try to check out an existing local branch called TEMPLATE
            try:
                self.repo.git.checkout("TEMPLATE")
            except GitCommandError:
                raise SyncExceptionError("Could not check out branch 'origin/TEMPLATE' or 'TEMPLATE'")

    def delete_tracked_template_branch_files(self):
        """
        Delete all tracked files and subsequent empty directories in the TEMPLATE branch
        """
        # Delete tracked files
        log.info("Deleting tracked files in 'TEMPLATE' branch")
        self._delete_tracked_files()
        self._clean_up_empty_dirs()

    def _delete_tracked_files(self):
        """
        Delete all tracked files in the repository
        """
        for the_file in self._get_tracked_files():
            file_path = Path(self.pipeline_dir) / the_file
            log.debug(f"Deleting {file_path}")
            try:
                file_path.unlink()
            except Exception as e:
                raise SyncExceptionError(e)

    def _clean_up_empty_dirs(self):
        """
        Delete empty directories in the repository

        Walks the directory tree from the bottom up, deleting empty directories as it goes.
        """
        # Track deleted child directories so we know they've been deleted when evaluating if the parent is empty
        deleted = set()

        for curr_dir, sub_dirs, files in os.walk(self.pipeline_dir, topdown=False):
            # Don't delete the root directory (should never happen due to .git, but just in case)
            if curr_dir == str(self.pipeline_dir):
                continue

            subdir_set = set(Path(curr_dir) / d for d in sub_dirs)
            currdir_is_empty = (len(subdir_set - deleted) == 0) and (len(files) == 0)
            if currdir_is_empty:
                log.debug(f"Deleting empty directory {curr_dir}")
                try:
                    Path(curr_dir).rmdir()
                except Exception as e:
                    raise SyncExceptionError(e)
                deleted.add(Path(curr_dir))

    def make_template_pipeline(self):
        """
        Delete all files and make a fresh template using the workflow variables
        """
        log.info("Making a new template pipeline using pipeline variables")

        # Only show error messages from pipeline creation
        logging.getLogger("nf_core.pipelines.create").setLevel(logging.ERROR)
        assert self.config_yml_path is not None
        assert self.config_yml is not None

        # Re-write the template yaml info from .nf-core.yml config
        if self.config_yml.template is not None:
            # Set force true in config to overwrite existing files

            self.config_yml.template.force = True
            with open(self.config_yml_path, "w") as config_path:
                yaml.safe_dump(self.config_yml.model_dump(exclude_none=True), config_path)

        try:
            pipeline_create_obj = nf_core.pipelines.create.create.PipelineCreate(
                outdir=str(self.pipeline_dir),
                from_config_file=True,
                no_git=True,
                force=True,
            )
            pipeline_create_obj.init_pipeline()

            # set force to false to avoid overwriting files in the future
            if self.config_yml.template is not None:
                self.config_yml.template = pipeline_create_obj.config
                # Set force true in config to overwrite existing files
                self.config_yml.template.force = False
                # Set outdir as the current directory to avoid local info leaking
                self.config_yml.template.outdir = "."
                # Update nf-core version
                self.config_yml.nf_core_version = nf_core.__version__
                dump_yaml_with_prettier(self.config_yml_path, self.config_yml.model_dump(exclude_none=True))

        except Exception as err:
            # Reset to where you were to prevent git getting messed up.
            self.repo.git.reset("--hard")
            raise SyncExceptionError(f"Failed to rebuild pipeline from template with error:\n{err}")

    def commit_template_changes(self):
        """If we have any changes with the new template files, make a git commit"""
        # Check that we have something to commit
        if not self.repo.is_dirty(untracked_files=True):
            log.info("Template contains no changes - no new commit created")
            return False
        # Commit changes
        try:
            newly_ignored_files = self._get_ignored_files()
            # add and commit all files except self.ignored_files
            # :! syntax to exclude files using git pathspec
            self.repo.git.add([f":!{f}" for f in self.ignored_files if f not in newly_ignored_files], all=True)
            self.repo.index.commit(f"Template update for nf-core/tools version {nf_core.__version__}")
            self.made_changes = True
            log.info("Committed changes to 'TEMPLATE' branch")
        except Exception as e:
            raise SyncExceptionError(f"Could not commit changes to TEMPLATE:\n{e}")
        return True

    def push_template_branch(self):
        """If we made any changes, push the TEMPLATE branch to the default remote
        and try to make a PR. If we don't have the auth token, try to figure out a URL
        for the PR and print this to the console.
        """
        log.info(f"Pushing TEMPLATE branch to remote: '{os.path.basename(self.pipeline_dir)}'")
        try:
            self.repo.git.push()
        except GitCommandError as e:
            raise PullRequestExceptionError(f"Could not push TEMPLATE branch:\n  {e}")

    def create_merge_base_branch(self):
        """Create a new branch from the updated TEMPLATE branch
        This branch will then be used to create the PR
        """
        # Check if branch exists already
        branch_list = [b.name for b in self.repo.branches]
        if self.merge_branch in branch_list:
            merge_branch_format = re.compile(rf"{self.original_merge_branch}-(\d+)")
            max_branch = max(
                [1]
                + [
                    int(merge_branch_format.match(branch).groups()[0])
                    for branch in branch_list
                    if merge_branch_format.match(branch)
                ]
            )
            new_branch = f"{self.original_merge_branch}-{max_branch + 1}"
            log.info(f"Branch already existed: '{self.merge_branch}', creating branch '{new_branch}' instead.")
            self.merge_branch = new_branch

        # Create new branch and checkout
        log.info(f"Checking out merge base branch '{self.merge_branch}'")
        try:
            self.repo.create_head(self.merge_branch)
        except GitCommandError as e:
            raise SyncExceptionError(f"Could not create new branch '{self.merge_branch}'\n{e}")

    def push_merge_branch(self):
        """Push the newly created merge branch to the remote repository"""
        log.info(f"Pushing '{self.merge_branch}' branch to remote")
        try:
            origin = self.repo.remote()
            origin.push(self.merge_branch)
        except GitCommandError as e:
            raise PullRequestExceptionError(f"Could not push branch '{self.merge_branch}':\n  {e}")

    def make_pull_request(self):
        """Create a pull request to a base branch (default: dev),
        from a head branch (default: TEMPLATE)

        Returns: An instance of class requests.Response
        """
        log.info("Submitting a pull request via the GitHub API")

        pr_title = f"Important! Template update for nf-core/tools v{nf_core.__version__}"
        blog_post_sentence = (
            f"For more details, check out the blog post: {self.blog_post}" if self.blog_post != "" else ""
        )
        pr_body_text = (
            "Version `{tag}` of [nf-core/tools](https://github.com/nf-core/tools) has just been released with updates to the nf-core template. "
            f"{blog_post_sentence}\n\n"
            "Please make sure to merge this pull-request as soon as possible, "
            f"resolving any merge conflicts in the `{self.merge_branch}` branch (or your own fork, if you prefer). "
            "Once complete, make a new minor release of your pipeline.\n\n"
            "For instructions on how to merge this PR, please see "
            "[https://nf-co.re/docs/contributing/sync/](https://nf-co.re/docs/contributing/sync/#merging-automated-prs).\n\n"
            "For more information about this release of [nf-core/tools](https://github.com/nf-core/tools), "
            "please see the `v{tag}` [release page](https://github.com/nf-core/tools/releases/tag/{tag})."
            "\n\n> [!NOTE]\n"
            "> Since nf-core/tools 3.5.0, older template update PRs will not be automatically closed, but will remain open in your pipeline repository. "
            "Older template PRs will be automatically closed once a newer template PR has been merged."
        ).format(tag=nf_core.__version__)

        # Make new pull-request
        stderr = rich.console.Console(stderr=True, force_terminal=nf_core.utils.rich_force_colors())
        with self.gh_api.cache_disabled():
            try:
                r = self.gh_api.request_retry(
                    f"https://api.github.com/repos/{self.gh_repo}/pulls",
                    post_data={
                        "title": pr_title,
                        "body": pr_body_text,
                        "maintainer_can_modify": True,
                        "head": self.merge_branch,
                        "base": self.from_branch,
                    },
                )
            except Exception as e:
                stderr.print_exception()
                raise PullRequestExceptionError(f"Something went badly wrong - {e}")
            else:
                self.gh_pr_returned_data = r.json()
                self.pr_url = self.gh_pr_returned_data["html_url"]
                log.debug(f"GitHub API PR worked, return code {r.status_code}")
                log.info(f"GitHub PR created: {self.gh_pr_returned_data['html_url']}")

    @staticmethod
    def _parse_json_response(response) -> tuple[Any, str]:
        """Helper method to parse JSON response and create pretty-printed string.

        Args:
            response: requests.Response object

        Returns:
            Tuple of (parsed_json, pretty_printed_str)
        """
        try:
            json_data = json.loads(response.content)
            return json_data, json.dumps(json_data, indent=4)
        except Exception:
            return response.content, str(response.content)

    def reset_target_dir(self):
        """
        Reset the target pipeline directory. Check out the original branch.
        """
        log.info(f"Checking out original branch: '{self.original_branch}'")
        try:
            self.repo.git.checkout(self.original_branch)
        except GitCommandError as e:
            raise SyncExceptionError(f"Could not reset to original branch `{self.original_branch}`:\n{e}")

    def _get_ignored_files(self) -> list[str]:
        """
        Get a list of all files in the repo ignored by git.
        """
        # -z separates with \0 and makes sure special characters are handled correctly
        raw_ignored_files = self.repo.git.ls_files(z=True, ignored=True, others=True, exclude_standard=True)
        return raw_ignored_files.split("\0")[:-1] if raw_ignored_files else []

    def _get_tracked_files(self) -> list[str]:
        """
        Get a list of all files in the repo tracked by git.
        """
        # -z separates with \0 and makes sure special characters are handled correctly
        raw_tracked_files = self.repo.git.ls_files(z=True)
        return raw_tracked_files.split("\0")[:-1] if raw_tracked_files else []
