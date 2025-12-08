from textual.app import ComposeResult
from textual.containers import Center
from textual.screen import Screen
from textual.widgets import Button, Footer, Header, Markdown, Static

from nf_core.utils import nfcore_logo

markdown = """
# Welcome to the nf-core pipeline creation wizard

This app will help you create a new Nextflow pipeline
from the [nf-core/tools pipeline template](https://github.com/nf-core/tools).

The template helps anyone benefit from nf-core best practices,
and is a requirement for nf-core pipelines.

> 💡 If you want to add a pipeline to nf-core, please
> [join on Slack](https://nf-co.re/join) and discuss your plans with the
> community as early as possible; _**ideally before you start on your pipeline!**_
> See the [nf-core guidelines](https://nf-co.re/docs/contributing/guidelines)
> and the [#new-pipelines](https://nfcore.slack.com/channels/new-pipelines)
> Slack channel for more information.
"""


class WelcomeScreen(Screen):
    """A welcome screen for the app."""

    def compose(self) -> ComposeResult:
        yield Header()
        yield Footer()
        yield Static(
            "\n" + "\n".join(nfcore_logo) + "\n",
            id="logo",
        )
        yield Markdown(markdown)
        yield Center(Button("Let's go!", id="start", variant="success"), classes="cta")
