from textwrap import dedent

from textual.app import ComposeResult
from textual.containers import Center
from textual.screen import Screen
from textual.widgets import Button, Footer, Header, Markdown, Static

from nf_core.pipelines.create.utils import add_hide_class
from nf_core.utils import nfcore_logo


class LoggingScreen(Screen):
    """A screen to show the final logs."""

    def compose(self) -> ComposeResult:
        yield Header()
        yield Footer()
        yield Markdown(
            dedent(
                """
                # Logging
                """
            )
        )
        yield Static(
            "\n" + "\n".join(nfcore_logo) + "\n",
            id="logo",
        )
        yield Markdown("Creating...")
        yield Center(self.parent.LOG_HANDLER.console)
        yield Center(
            Button("Back", id="back", variant="default", classes="hide"),
            Button("Continue", id="close_screen", variant="success", classes="hide"),
            Button("Continue", id="exit", variant="success", classes="hide"),
            Button("Close App", id="close_app", variant="success", classes="hide"),
            classes="cta",
        )

    def on_screen_resume(self):
        """Hide all buttons as disabled on screen resume."""
        button_ids = ["back", "close_screen", "exit", "close_app"]
        for button in self.query("Button"):
            if button.id in button_ids:
                add_hide_class(self.parent, button.id)

    def on_screen_suspend(self):
        """Clear console on screen suspend."""
        self.parent.LOG_HANDLER.console.clear()
