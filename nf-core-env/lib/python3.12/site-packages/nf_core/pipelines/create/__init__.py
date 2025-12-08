"""A Textual app to create a pipeline."""

import logging

import click
from rich.logging import RichHandler
from textual.app import App
from textual.widgets import Button, Switch

from nf_core.pipelines.create import utils
from nf_core.pipelines.create.basicdetails import BasicDetails
from nf_core.pipelines.create.custompipeline import CustomPipeline
from nf_core.pipelines.create.finaldetails import FinalDetails
from nf_core.pipelines.create.githubexit import GithubExit
from nf_core.pipelines.create.githubrepo import GithubRepo
from nf_core.pipelines.create.githubrepoquestion import GithubRepoQuestion
from nf_core.pipelines.create.loggingscreen import LoggingScreen
from nf_core.pipelines.create.nfcorepipeline import NfcorePipeline
from nf_core.pipelines.create.pipelinetype import ChoosePipelineType
from nf_core.pipelines.create.welcome import WelcomeScreen

logger = logging.getLogger(__name__)
rich_log_handler = RichHandler(
    console=utils.LoggingConsole(classes="log_console"),
    level=logging.INFO,
    rich_tracebacks=True,
    show_time=False,
    show_path=False,
    markup=True,
    tracebacks_suppress=[click],
)
logger.addHandler(rich_log_handler)


class PipelineCreateApp(App[utils.CreateConfig]):
    """A Textual app to manage stopwatches."""

    CSS_PATH = "create.tcss"
    TITLE = "nf-core pipelines create"
    SUB_TITLE = "Create a new pipeline with the nf-core pipeline template"
    BINDINGS = [
        ("d", "toggle_dark", "Toggle dark mode"),
        ("q", "quit", "Quit"),
        ("a", "toggle_all", "Toggle all"),
    ]
    SCREENS = {
        "welcome": WelcomeScreen,
        "basic_details": BasicDetails,
        "choose_type": ChoosePipelineType,
        "type_custom": CustomPipeline,
        "type_nfcore": NfcorePipeline,
        "final_details": FinalDetails,
        "logging": LoggingScreen,
        "github_repo_question": GithubRepoQuestion,
        "github_repo": GithubRepo,
        "github_exit": GithubExit,
    }

    # Initialise config as empty
    TEMPLATE_CONFIG = utils.CreateConfig()

    # Initialise pipeline type
    NFCORE_PIPELINE = True

    # Log handler
    LOG_HANDLER = rich_log_handler
    # Logging state
    LOGGING_STATE = None

    # Template features
    template_features_yml = utils.load_features_yaml()

    def on_mount(self) -> None:
        self.push_screen("welcome")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle all button pressed events."""
        if event.button.id == "start":
            self.push_screen("choose_type")
        elif event.button.id == "type_nfcore":
            self.NFCORE_PIPELINE = True
            utils.NFCORE_PIPELINE_GLOBAL = True
            self.push_screen("basic_details")
        elif event.button.id == "type_custom":
            self.NFCORE_PIPELINE = False
            utils.NFCORE_PIPELINE_GLOBAL = False
            self.push_screen("basic_details")
        elif event.button.id == "continue":
            self.push_screen("final_details")
        elif event.button.id == "github_repo":
            self.push_screen("github_repo")
        elif event.button.id == "close_screen":
            self.push_screen("github_repo_question")
        elif event.button.id == "exit":
            self.push_screen("github_exit")
        if event.button.id == "close_app":
            self.exit(return_code=0)
        if event.button.id == "back":
            self.pop_screen()

    def action_toggle_dark(self) -> None:
        """An action to toggle dark mode."""
        self.theme: str = "textual-dark" if self.theme == "textual-light" else "textual-light"

    def action_toggle_all(self) -> None:
        """An action to toggle all Switches."""
        switches = self.screen.query("Switch")
        if not switches:
            return  # No Switches widgets found
        # Determine the new state based on the first switch
        new_state = not switches.first().value if switches.first() else True
        for switch in switches:
            switch.value = new_state
        self.refresh()
