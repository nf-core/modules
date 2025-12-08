"""
The ComponentsTest class handles the generation and testing of nf-test snapshots.
"""

import logging
import os
import re
from pathlib import Path

import questionary
from rich import print
from rich.panel import Panel
from rich.prompt import Confirm
from rich.syntax import Syntax
from rich.text import Text

import nf_core.utils
from nf_core.components.components_command import ComponentCommand

log = logging.getLogger(__name__)


class ComponentsTest(ComponentCommand):  # type: ignore[misc]
    """
    Class to generate and test nf-test snapshots for modules.

    ...

    Attributes
    ----------
    component_type : str
        type of component to test (modules or subworkflows)
    component_name : str
        name of the tool to run tests for
    directory: str
        path to modules repository directory
    no_prompts : bool
        flat indicating if prompts are used
    remote_url : str
        URL of the remote repository
    branch : str
        branch of the remote repository
    verbose : bool
        flag indicating if verbose output should be used
    update : bool
        flag indicating if the existing snapshot should be updated
    once : bool
        flag indicating if the test should be run only once
    profile : str
        container software to use (docker, singularity or conda)

    Methods
    -------
    run():
        Run test steps
    check_inputs():
        Check inputs. Ask for component_name if not provided and check that the directory exists
    generate_snapshot():
        Generate the nf-test snapshot using `nf-test test` command
    check_snapshot_stability():
        Run the nf-test twice and check if the snapshot changes
    """

    def __init__(
        self,
        component_type: str,
        component_name: str | None = None,
        directory: str = ".",
        no_prompts: bool = False,
        remote_url: str | None = None,
        branch: str | None = None,
        verbose: bool = False,
        update: bool = False,
        once: bool = False,
        profile: str | None = None,
    ):
        super().__init__(component_type, directory, remote_url, branch, no_prompts=no_prompts)
        self.component_name = component_name
        self.remote_url = remote_url
        self.branch = branch
        self.errors: list[str] = []
        self.verbose = verbose
        self.obsolete_snapshots: bool = False
        self.update = update
        self.once = once
        self.profile = profile

    def run(self) -> None:
        """Run build steps"""
        self.check_inputs()
        os.environ["NFT_DIFF"] = "pdiff"  # set nf-test differ to pdiff to get a better diff output
        os.environ["NFT_DIFF_ARGS"] = (
            "--line-numbers --expand-tabs=2"  # taken from https://code.askimed.com/nf-test/docs/assertions/snapshots/#snapshot-differences
        )
        with nf_core.utils.set_wd(self.directory):
            self.check_snapshot_stability()
        if len(self.errors) > 0:
            errors = "\n - ".join(self.errors)
            raise UserWarning(f"Ran, but found errors:\n - {errors}")
        else:
            log.info("All tests passed!")

    def check_inputs(self) -> None:
        """Do more complex checks about supplied flags."""
        # Check modules directory structure
        if self.component_type == "modules":
            self.check_modules_structure()

        # Get the component name if not specified
        if self.component_name is None:
            if self.no_prompts:
                raise UserWarning(
                    f"{self.component_type[:-1].title()} name not provided and prompts deactivated. Please provide the {self.component_type[:-1]} name{' as TOOL/SUBTOOL or TOOL' if self.component_type == 'modules' else ''}."
                )
            else:
                try:
                    self.component_name = questionary.autocomplete(
                        "Tool name:" if self.component_type == "modules" else "Subworkflow name:",
                        choices=self.components_from_repo(self.org),
                        style=nf_core.utils.nfcore_question_style,
                    ).unsafe_ask()
                except LookupError:
                    raise

        self.component_dir = Path(self.component_type, self.modules_repo.repo_path, *self.component_name.split("/"))

        # First, sanity check that the module directory exists
        if not Path(self.directory, self.component_dir).is_dir():
            raise UserWarning(
                f"Cannot find directory '{self.component_dir}'.{' Should be TOOL/SUBTOOL or TOOL' if self.component_type == 'modules' else ''}"
            )

        # Check container software to use
        if os.environ.get("PROFILE") is None and self.profile is None:
            os.environ["PROFILE"] = ""
            if self.no_prompts:
                log.info(
                    "Setting environment variable '$PROFILE' to Docker as not set otherwise.\n"
                    "To use Singularity set 'export PROFILE=singularity' in your shell before running this command."
                )
                os.environ["PROFILE"] = "docker"
            else:
                question = {
                    "type": "list",
                    "name": "profile",
                    "message": "Choose container software to run the test with",
                    "choices": ["Docker", "Singularity", "Conda"],
                }
                answer = questionary.unsafe_prompt([question], style=nf_core.utils.nfcore_question_style)
                profile = answer["profile"].lower()
                os.environ["PROFILE"] = profile

    def display_nftest_output(self, nftest_out: bytes, nftest_err: bytes) -> None:
        nftest_output = Text.from_ansi(nftest_out.decode())
        print(Panel(nftest_output, title="nf-test output"))
        if nftest_err:
            syntax = Syntax(nftest_err.decode(), "diff", theme="ansi_dark")
            print(Panel(syntax, title="nf-test error"))
            if "Different Snapshot:" in nftest_err.decode():
                log.error("nf-test failed due to differences in the snapshots")
                # prompt to update snapshot
                if self.no_prompts:
                    log.info("Updating snapshot")
                    self.update = True
                elif self.update is None:
                    answer = Confirm.ask(
                        "[bold][blue]?[/] nf-test found differences in the snapshot. Do you want to update it?",
                        default=True,
                    )
                    if answer:
                        log.info("Updating snapshot")
                        self.update = True
                    else:
                        log.debug("Snapshot not updated")
                if self.update:
                    # update snapshot using nf-test --update-snapshot
                    self.generate_snapshot()

            else:
                self.errors.append("nf-test failed")

    def generate_snapshot(self) -> bool:
        """Generate the nf-test snapshot using `nf-test test` command

        returns True if the test was successful, False otherwise
        """

        log.debug("Running nf-test test")

        # set verbose flag if self.verbose is True
        verbose = "--verbose --debug" if self.verbose else ""
        update = "--update-snapshot" if self.update else ""
        self.update = False  # reset self.update to False to test if the new snapshot is stable
        tag = f"subworkflows/{self.component_name}" if self.component_type == "subworkflows" else self.component_name
        profile = self.profile if self.profile else os.environ["PROFILE"]

        result = nf_core.utils.run_cmd(
            "nf-test",
            f"test --tag {tag} --profile {profile} {verbose} {update}",
        )
        if result is not None:
            nftest_out, nftest_err = result
            self.display_nftest_output(nftest_out, nftest_err)
            # check if nftest_out contains obsolete snapshots
            pattern = r"Snapshot Summary:.*?(\d+)\s+obsolete"
            compiled_pattern = re.compile(pattern, re.DOTALL)  # re.DOTALL to allow . to match newlines
            obsolete_snapshots = compiled_pattern.search(nftest_out.decode())
            if obsolete_snapshots:
                self.obsolete_snapshots = True
            # check if nf-test was successful
            if "Assertion failed:" in nftest_out.decode():
                return False
            elif "No tests to execute." in nftest_out.decode():
                log.error("Nothing to execute. Is the file 'main.nf.test' missing?")
                self.errors.append("Nothing to execute. Is the file 'main.nf.test' missing?")
                return False
            else:
                log.debug("nf-test successful")
                return True
        else:
            log.error("nf-test failed")
            self.errors.append("nf-test failed")
            return False

    def check_snapshot_stability(self) -> bool:
        """Run the nf-test twice and check if the snapshot changes"""
        log.info("Generating nf-test snapshot")
        if not self.generate_snapshot():
            return False  # stop here if the first run failed
        elif self.once:
            return True  # stop here if the test should be run only once
        log.info("Generating nf-test snapshot again to check stability")
        if not self.generate_snapshot():
            log.error("nf-test snapshot is not stable")
            self.errors.append("nf-test snapshot is not stable")
            return False

        else:
            if self.obsolete_snapshots:
                # ask if the user wants to remove obsolete snapshots using nf-test --clean-snapshot
                if self.no_prompts or Confirm.ask(
                    "nf-test found obsolete snapshots. Do you want to remove them?", default=True
                ):
                    profile = self.profile if self.profile else os.environ["PROFILE"]
                    log.info("Removing obsolete snapshots")
                    nf_core.utils.run_cmd(
                        "nf-test",
                        f"test --tag {self.component_name} --profile {profile} --clean-snapshot",
                    )
                else:
                    log.debug("Obsolete snapshots not removed")
            return True
