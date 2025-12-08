"""A Textual app to create a pipeline."""

from pathlib import Path
from textwrap import dedent

from textual import on, work
from textual.app import ComposeResult
from textual.containers import Center, Horizontal
from textual.screen import Screen
from textual.widgets import Button, Footer, Header, Input, Markdown

from nf_core.pipelines.create.create import PipelineCreate
from nf_core.pipelines.create.utils import ShowLogs, TextInput, add_hide_class, remove_hide_class

pipeline_exists_warn = """
> ⚠️  **The pipeline you are trying to create already exists.**
>
> If you continue, you will **override** the existing pipeline.
> Please change the pipeline or organisation name to create a different pipeline.
> Alternatively, provide a different output directory.
"""


class FinalDetails(Screen):
    """Name, description, author, etc."""

    def compose(self) -> ComposeResult:
        yield Header()
        yield Footer()
        yield Markdown(
            dedent(
                """
                # Final details
                """
            )
        )

        with Horizontal():
            yield TextInput(
                "version",
                "Version",
                "First version of the pipeline",
                "1.0.0dev",
                classes="column",
            )
            yield TextInput(
                "outdir",
                "Output directory",
                "Path to the output directory where the pipeline will be created",
                ".",
                classes="column",
            )

        yield Markdown(dedent(pipeline_exists_warn), id="exist_warn", classes="hide")

        yield Center(
            Button("Back", id="back", variant="default"),
            Button("Finish", id="finish", variant="success"),
            classes="cta",
        )

    @on(Button.Pressed, "#finish")
    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Save fields to the config."""
        new_config = {}
        for text_input in self.query("TextInput"):
            this_input = text_input.query_one(Input)
            validation_result = this_input.validate(this_input.value)
            new_config[text_input.field_id] = this_input.value
            if not validation_result.is_valid:
                text_input.query_one(".validation_msg").update("\n".join(validation_result.failure_descriptions))
            else:
                text_input.query_one(".validation_msg").update("")
        try:
            self.parent.TEMPLATE_CONFIG.__dict__.update(new_config)
        except ValueError:
            pass

        # Create the new pipeline
        self._create_pipeline()
        self.parent.LOGGING_STATE = "pipeline created"
        self.parent.push_screen("logging")

    @on(Input.Changed)
    @on(Input.Submitted)
    def show_exists_warn(self):
        """Check if the pipeline exists on every input change or submitted.
        If the pipeline exists, show warning message saying that it will be overridden."""
        outdir = ""
        for text_input in self.query("TextInput"):
            this_input = text_input.query_one(Input)
            if text_input.field_id == "outdir":
                outdir = this_input.value
        if Path(outdir, self.parent.TEMPLATE_CONFIG.org + "-" + self.parent.TEMPLATE_CONFIG.name).is_dir():
            remove_hide_class(self.parent, "exist_warn")

    def on_screen_resume(self):
        """Hide warn message on screen resume."""
        add_hide_class(self.parent, "exist_warn")

    @work(thread=True, exclusive=True)
    def _create_pipeline(self) -> None:
        """Create the pipeline."""
        self.post_message(ShowLogs())
        create_obj = PipelineCreate(
            template_config=self.parent.TEMPLATE_CONFIG,
            is_interactive=True,
        )
        create_obj.init_pipeline()
        remove_hide_class(self.parent, "close_screen")
