from textwrap import dedent

from textual.app import ComposeResult
from textual.containers import Center
from textual.screen import Screen
from textual.widgets import Button, Footer, Header, Markdown, Static

from nf_core.utils import nfcore_logo

exit_help_text_markdown = """
If you would like to create the GitHub repository later, you can do it manually by following these steps:

1. Create a new GitHub repository
2. Add the remote to your local repository:
    ```bash
    cd <pipeline_directory>
    git remote add origin git@github.com:<username>/<repo_name>.git
    ```
3. Push the code to the remote:
    ```bash
    git push --all origin
    ```
    > 💡 Note the `--all` flag: this is needed to push all branches to the remote.
"""


class GithubExit(Screen):
    """A screen to show a help text when a GitHub repo is NOT created."""

    def compose(self) -> ComposeResult:
        yield Header()
        yield Footer()
        yield Markdown(
            dedent(
                """
                # HowTo create a GitHub repository
                """
            )
        )
        yield Static(
            "\n" + "\n".join(nfcore_logo) + "\n",
            id="logo",
        )
        yield Markdown(exit_help_text_markdown)
        yield Center(
            Button("Close", id="close_app", variant="success"),
            classes="cta",
        )
