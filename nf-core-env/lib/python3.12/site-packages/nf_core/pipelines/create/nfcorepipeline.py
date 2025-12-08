from textwrap import dedent

from textual import on
from textual.app import ComposeResult
from textual.containers import Center, ScrollableContainer
from textual.screen import Screen
from textual.widgets import Button, Footer, Header, Markdown, Switch

from nf_core.pipelines.create.utils import PipelineFeature


class NfcorePipeline(Screen):
    """Select if the pipeline will use genomic data."""

    def compose(self) -> ComposeResult:
        yield Header()
        yield Footer()
        yield Markdown(
            dedent(
                """
                # Template features
                """
            )
        )
        yield ScrollableContainer(id="features")
        yield Center(
            Button("Back", id="back", variant="default"),
            Button("Continue", id="continue", variant="success"),
            classes="cta",
        )

    def on_mount(self) -> None:
        for section_name, section in self.parent.template_features_yml.items():
            section_title = section["name"]
            features = section["features"]
            show_section = False
            self.query_one("#features").mount(
                Markdown(section_title, id=section_name),
            )
            for name, feature in features.items():
                if feature["nfcore_pipelines"]:
                    show_section = True
                    self.query_one("#features").mount(
                        PipelineFeature(
                            feature["help_text"],
                            feature["short_description"],
                            feature["description"],
                            name,
                            feature["default"],
                        )
                    )
            if not show_section:
                self.query_one("#features").query_one(f"#{section_name}").remove()

    @on(Button.Pressed, "#continue")
    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Save fields to the config."""
        skip = []
        for feature_input in self.query("PipelineFeature"):
            this_switch = feature_input.query_one(Switch)
            if not this_switch.value:
                skip.append(this_switch.id)
        self.parent.TEMPLATE_CONFIG.__dict__.update({"skip_features": skip, "is_nfcore": True})
