from textual.app import ComposeResult
from textual.containers import Center, Grid
from textual.screen import Screen
from textual.widgets import Button, Footer, Header, Markdown

markdown_intro = """
# Choose pipeline type
"""

markdown_type_nfcore = """
## Choose _"nf-core"_ if:

* You want your pipeline to be part of the nf-core community
* Your pipeline has been accepted via the [nf-core proposal process](https://github.com/nf-core/proposals)
"""
markdown_type_custom = """
## Choose _"Custom"_ if:

* Your pipeline will _never_ be part of nf-core
* You want full control over *all* features that are included from the template
      (including those that are mandatory for nf-core).
"""

markdown_details = """
## What's the difference?

Choosing _"nf-core"_ effectively pre-selects the following template features:

* GitHub Actions continuous-integration configuration files:
    * Pipeline test runs: Small-scale (GitHub) and large-scale (AWS)
    * Code formatting checks with [Prettier](https://prettier.io/)
    * Auto-fix linting functionality using [@nf-core-bot](https://github.com/nf-core-bot)
    * Marking old issues as stale
* Inclusion of [shared nf-core configuration profiles](https://nf-co.re/configs)
"""


class ChoosePipelineType(Screen):
    """Choose whether this will be an nf-core pipeline or not."""

    def compose(self) -> ComposeResult:
        yield Header()
        yield Footer()
        yield Markdown(markdown_intro)
        yield Grid(
            Center(
                Markdown(markdown_type_nfcore),
                Center(Button("nf-core", id="type_nfcore", variant="success")),
            ),
            Center(
                Markdown(markdown_type_custom),
                Center(Button("Custom", id="type_custom", variant="primary")),
            ),
            classes="col-2 pipeline-type-grid",
        )
        yield Markdown(markdown_details)
