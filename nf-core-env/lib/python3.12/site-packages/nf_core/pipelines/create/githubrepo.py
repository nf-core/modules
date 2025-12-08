import logging
import os
from pathlib import Path
from textwrap import dedent

import git
import yaml
from github import Github, GithubException, UnknownObjectException
from rich.text import Text
from textual import on, work
from textual.app import ComposeResult
from textual.containers import Center, Horizontal, Vertical
from textual.message import Message
from textual.screen import Screen
from textual.widgets import Button, Footer, Header, Input, Markdown, Static, Switch

from nf_core.pipelines.create.utils import ShowLogs, TextInput, remove_hide_class

log = logging.getLogger(__name__)

github_org_help = """
> ⚠️  **You can't create a repository directly in the nf-core organisation.**
>
> Please create the pipeline repo to an organisation where you have access or use your user account.
> A core-team member will be able to transfer the repo to nf-core once the development has started.

> 💡 Your GitHub user account will be used by default if `nf-core` is given as the org name.
"""


class GithubRepo(Screen):
    """Create a GitHub repository and push all branches."""

    def compose(self) -> ComposeResult:
        yield Header()
        yield Footer()
        gh_user, gh_token = self._get_github_credentials()
        github_text_markdown = dedent(
            """
            # Create GitHub repository

            Now that we have created a new pipeline locally, we can
            create a new GitHub repository and push the code to it.
            """
        )
        if gh_user:
            github_text_markdown += f">\n> 💡 _Found GitHub username {'and token ' if gh_token else ''}in local [GitHub CLI](https://cli.github.com/) config_\n>\n"
        yield Markdown(github_text_markdown)
        with Horizontal(classes="ghrepo-cols"):
            yield TextInput(
                "gh_username",
                "GitHub username",
                "Your GitHub username",
                default=gh_user[0] if gh_user is not None else "GitHub username",
                classes="column",
            )
            yield TextInput(
                "token",
                "GitHub token",
                Text.from_markup("Your GitHub [link=https://shorturl.at/RKJzS]personal access token[/link] for login."),
                default=gh_token if gh_token is not None else "GitHub token",
                password=True,
                classes="column",
            )
            yield Button("Show", id="show_password")
            yield Button("Hide", id="hide_password")
        with Horizontal(classes="ghrepo-cols"):
            yield TextInput(
                "repo_org",
                "Organisation name",
                "The name of the organisation where the GitHub repo will be created",
                default=self.parent.TEMPLATE_CONFIG.org,
                classes="column",
            )
            yield TextInput(
                "repo_name",
                "Repository name",
                "The name of the new GitHub repository",
                default=self.parent.TEMPLATE_CONFIG.name,
                classes="column",
            )
        if self.parent.TEMPLATE_CONFIG.is_nfcore:
            yield Markdown(dedent(github_org_help))
        with Horizontal(classes="ghrepo-cols"):
            yield Switch(value=False, id="private")
            with Vertical():
                yield Static("Private", classes="")
                yield Static("Select to make the new GitHub repo private.", classes="feature_subtitle")
        yield Center(
            Button("Back", id="back", variant="default"),
            Button("Create GitHub repo", id="create_github", variant="success"),
            Button("Finish without creating a repo", id="exit", variant="primary"),
            classes="cta",
        )

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Create a GitHub repo or show help message and exit"""
        if event.button.id == "show_password":
            self.add_class("displayed")
            text_input = self.query_one("#token", TextInput)
            text_input.query_one(Input).password = False
        elif event.button.id == "hide_password":
            self.remove_class("displayed")
            text_input = self.query_one("#token", TextInput)
            text_input.query_one(Input).password = True
        elif event.button.id == "create_github":
            # Create a GitHub repo

            # Save GitHub username, token and repo name
            github_variables = {}
            for text_input in self.query("TextInput"):
                this_input = text_input.query_one(Input)
                github_variables[text_input.field_id] = this_input.value
            # Save GitHub repo config
            for switch_input in self.query("Switch"):
                github_variables[switch_input.id] = switch_input.value

            # Pipeline git repo
            pipeline_repo = git.Repo.init(
                Path(self.parent.TEMPLATE_CONFIG.outdir)
                / Path(self.parent.TEMPLATE_CONFIG.org + "-" + self.parent.TEMPLATE_CONFIG.name)
            )

            # GitHub authentication
            if github_variables["token"]:
                github_auth = self._github_authentication(github_variables["gh_username"], github_variables["token"])
            else:
                raise UserWarning(
                    f"Could not authenticate to GitHub with user name '{github_variables['gh_username']}'."
                    "Please provide an authentication token or set the environment variable 'GITHUB_AUTH_TOKEN'."
                )

            user = github_auth.get_user()
            org = None
            # Make sure that the authentication was successful
            try:
                user.login
                log.debug("GitHub authentication successful")
            except GithubException:
                raise UserWarning(
                    f"Could not authenticate to GitHub with user name '{github_variables['gh_username']}'."
                    "Please make sure that the provided user name and token are correct."
                )

            # Check if organisation exists
            # If the organisation is nf-core or it doesn't exist, the repo will be created in the user account
            if github_variables["repo_org"] != "nf-core":
                try:
                    org = github_auth.get_organization(github_variables["repo_org"])
                    log.info(
                        f"Repo will be created in the GitHub organisation account '{github_variables['repo_org']}'"
                    )
                except UnknownObjectException:
                    log.warn(f"Provided organisation '{github_variables['repo_org']}' not found. ")

            # Create the repo
            try:
                if org:
                    self._create_repo_and_push(
                        org,
                        github_variables["repo_name"],
                        pipeline_repo,
                        github_variables["private"],
                    )
                else:
                    # Create the repo in the user's account
                    log.info(
                        f"Repo will be created in the GitHub organisation account '{github_variables['gh_username']}'"
                    )
                    self._create_repo_and_push(
                        user,
                        github_variables["repo_name"],
                        pipeline_repo,
                        github_variables["private"],
                    )
            except UserWarning as e:
                log.error(f"There was an error with message: {e}")
                self.parent.push_screen("github_exit")

            self.parent.LOGGING_STATE = "repo created"
            self.parent.push_screen("logging")

    class RepoExists(Message):
        """Custom message to indicate that the GitHub repo already exists."""

        pass

    @on(RepoExists)
    def show_github_info_button(self) -> None:
        remove_hide_class(self.parent, "exit")
        remove_hide_class(self.parent, "back")

    @work(thread=True, exclusive=True)
    def _create_repo_and_push(self, org, repo_name, pipeline_repo, private):
        """Create a GitHub repository and push all branches."""
        self.post_message(ShowLogs())
        # Check if repo already exists
        try:
            repo = org.get_repo(repo_name)
            # Check if it has a commit history
            try:
                repo.get_commits().totalCount
                raise UserWarning(f"GitHub repository '{repo_name}' already exists")
            except GithubException:
                # Repo is empty
                repo_exists = True
            except UserWarning as e:
                # Repo already exists
                log.error(e)
                self.post_message(self.RepoExists())
                return
        except UnknownObjectException:
            # Repo doesn't exist
            repo_exists = False

        # Create the repo
        if not repo_exists:
            repo = org.create_repo(repo_name, description=self.parent.TEMPLATE_CONFIG.description, private=private)
            log.info(f"GitHub repository '{repo_name}' created successfully")
            remove_hide_class(self.parent, "close_app")

        # Add the remote
        try:
            pipeline_repo.create_remote("origin", repo.clone_url)
        except git.exc.GitCommandError:
            # Remote already exists
            pass
        # Push all branches
        pipeline_repo.remotes.origin.push(all=True).raise_if_error()

    def _github_authentication(self, gh_username, gh_token):
        """Authenticate to GitHub"""
        log.debug(f"Authenticating GitHub as {gh_username}")
        github_auth = Github(gh_username, gh_token)
        return github_auth

    def _get_github_credentials(self):
        """Get GitHub credentials"""
        gh_user = None
        gh_token = None
        # Use gh CLI config if installed
        gh_cli_config_fn = os.path.expanduser("~/.config/gh/hosts.yml")
        if os.path.exists(gh_cli_config_fn):
            try:
                with open(gh_cli_config_fn) as fh:
                    gh_cli_config = yaml.safe_load(fh)
                    gh_user = (gh_cli_config["github.com"]["user"],)
                    gh_token = gh_cli_config["github.com"]["oauth_token"]
            except KeyError:
                pass
        # If gh CLI not installed, try to get credentials from environment variables
        elif os.environ.get("GITHUB_TOKEN") is not None:
            gh_token = self.auth = os.environ["GITHUB_TOKEN"]
        return (gh_user, gh_token)
