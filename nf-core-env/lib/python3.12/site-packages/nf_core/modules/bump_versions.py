"""
Bump versions for all modules on nf-core/modules
or for a single module
"""

import logging
import os
import re
from pathlib import Path

import questionary
import yaml
from rich.box import ROUNDED
from rich.console import Console
from rich.markdown import Markdown
from rich.panel import Panel
from rich.progress import BarColumn, Progress
from rich.table import Table

import nf_core.modules.modules_utils
import nf_core.utils
from nf_core.components.components_command import ComponentCommand
from nf_core.components.nfcore_component import NFCoreComponent
from nf_core.utils import NFCoreYamlConfig, custom_yaml_dumper, rich_force_colors
from nf_core.utils import plural_s as _s

log = logging.getLogger(__name__)


class ModuleVersionBumper(ComponentCommand):
    def __init__(
        self,
        pipeline_dir: str | Path,
        remote_url: str | None = None,
        branch: str | None = None,
        no_pull: bool = False,
    ):
        super().__init__("modules", pipeline_dir, remote_url, branch, no_pull)

        self.up_to_date: list[tuple[str, str]] = []
        self.updated: list[tuple[str, str]] = []
        self.failed: list[tuple[str, str]] = []
        self.ignored: list[tuple[str, str]] = []
        self.show_up_to_date: bool | None = None
        self.tools_config: NFCoreYamlConfig | None

    def bump_versions(
        self,
        module: str | None = None,
        all_modules: bool = False,
        show_up_to_date: bool = False,
        dry_run: bool = False,
    ) -> list[NFCoreComponent]:
        """
        Bump the container and conda version of single module or all modules.

        If module is the name of a directory in the modules directory, all modules in that directory will be bumped.

        Looks for a bioconda tool version in the `main.nf` file of the module and checks whether
        are more recent version is available. If yes, then tries to get docker/singularity
        container links and replace the bioconda version and the container links in the main.nf file
        of the respective module.

        Args:
            module: a specific module to update
            all_modules: whether to bump versions for all modules
            show_up_to_date: whether to show up-to-date modules as well
            dry_run: whether to dry run the command

        Returns:
            list[NFCoreComponent]: the updated modules
        """
        self.up_to_date = []
        self.updated = []
        self.failed = []
        self.ignored = []
        self.show_up_to_date = show_up_to_date

        # Check modules directory structure
        self.check_modules_structure()

        # Verify that this is not a pipeline
        if not self.repo_type == "modules":
            raise nf_core.modules.modules_utils.ModuleExceptionError(
                "This command only works on the nf-core/modules repository, not on pipelines!"
            )

        # Get list of all modules
        _, nfcore_modules = nf_core.modules.modules_utils.get_installed_modules(self.directory)

        # Load the .nf-core.yml config
        _, self.tools_config = nf_core.utils.load_tools_config(self.directory)

        # Prompt for module or all
        if module is None and not all_modules:
            question = {
                "type": "list",
                "name": "all_modules",
                "message": "Bump versions for all modules or a single named module?",
                "choices": ["All modules", "Named module"],
            }
            answer = questionary.unsafe_prompt([question], style=nf_core.utils.nfcore_question_style)
            if answer["all_modules"] == "All modules":
                all_modules = True
            else:
                module = questionary.autocomplete(
                    "Tool name:",
                    choices=[m.component_name for m in nfcore_modules],
                    style=nf_core.utils.nfcore_question_style,
                ).unsafe_ask()

        if module:
            self.show_up_to_date = True
            if all_modules:
                raise nf_core.modules.modules_utils.ModuleExceptionError(
                    "You cannot specify a tool and request all tools to be bumped."
                )
            nfcore_modules = nf_core.modules.modules_utils.filter_modules_by_name(nfcore_modules, module)

            if len(nfcore_modules) == 0:
                raise nf_core.modules.modules_utils.ModuleExceptionError(
                    f"Could not find the specified module: '{module}'"
                )

        # mainly used for testing, return the list of nfcore_modules selected
        if dry_run:
            return nfcore_modules

        progress_bar = Progress(
            "[bold blue]{task.description}",
            BarColumn(bar_width=None),
            "[magenta]{task.completed} of {task.total}[reset] » [bold yellow]{task.fields[test_name]}",
            transient=True,
            disable=os.environ.get("HIDE_PROGRESS", None) is not None,
        )
        with progress_bar:
            bump_progress = progress_bar.add_task(
                "Bumping nf-core modules versions",
                total=len(nfcore_modules),
                test_name=nfcore_modules[0].component_name,
            )
            for mod in nfcore_modules:
                progress_bar.update(bump_progress, advance=1, test_name=mod.component_name)
                self.bump_module_version(mod)

        self._print_results()

        return nfcore_modules

    def bump_module_version(self, module: NFCoreComponent) -> bool:
        """
        Bump the bioconda and container version of a single NFCoreComponent

        Args:
            module: NFCoreComponent
        """
        config_version = None
        bioconda_packages = []
        try:
            # Extract bioconda version from `environment.yml`
            bioconda_packages = self.get_bioconda_version(module)
        except FileNotFoundError:
            # try it in the main.nf instead
            try:
                with open(module.main_nf) as fh:
                    for line in fh:
                        if "bioconda::" in line:
                            bioconda_packages = [b for b in line.split() if "bioconda::" in b]
            except FileNotFoundError:
                log.error(
                    f"Neither `environment.yml` nor `main.nf` of {module.component_name} module could be read to get bioconada version of used tools."
                )

        # If multiple versions - don't update! (can't update mulled containers)
        if not bioconda_packages or len(bioconda_packages) > 1:
            self.failed.append(("Ignoring mulled container", module.component_name))
            return False

        # Don't update if blocked in blacklist
        bump_versions_config: dict[str, str] = getattr(self.tools_config, "bump-versions", {}) or {}
        if module.component_name in bump_versions_config:
            config_version = bump_versions_config[module.component_name]
            if not config_version:
                self.ignored.append(("Omitting module due to config.", module.component_name))
                return False

        # check for correct version and newer versions
        bioconda_tool_name = bioconda_packages[0].split("=")[0].replace("bioconda::", "").strip("'").strip('"')
        bp = bioconda_packages[0]
        bp = bp.strip("'").strip('"')
        bioconda_version = bp.split("=")[1]

        if not config_version:
            try:
                response = nf_core.utils.anaconda_package(bp)
            except (LookupError, ValueError):
                self.failed.append(
                    (
                        f"Conda version not specified correctly: {Path(module.main_nf).relative_to(self.directory)}",
                        module.component_name,
                    )
                )
                return False

            # Check that required version is available at all
            if bioconda_version not in response.get("versions"):
                self.failed.append((f"Conda package had unknown version: `{module.main_nf}`", module.component_name))
                return False

            # Check version is latest available
            last_ver = response.get("latest_version")
        else:
            last_ver = config_version

        if last_ver is not None and last_ver != bioconda_version:
            log.debug(f"Updating version for {module.component_name}")
            # Get docker and singularity container links
            try:
                docker_img, singularity_img = nf_core.utils.get_biocontainer_tag(bioconda_tool_name, last_ver)
            except LookupError as e:
                self.failed.append((f"Could not download container tags: {e}", module.component_name))
                return False

            patterns = [
                (rf"biocontainers/{bioconda_tool_name}:[^'\"\s]+", docker_img),
                (
                    rf"https://depot.galaxyproject.org/singularity/{bioconda_tool_name}:[^'\"\s]+",
                    singularity_img,
                ),
            ]

            with open(module.main_nf) as fh:
                content = fh.read()

            # Go over file content of main.nf and find replacements
            for pattern in patterns:
                found_match = False
                newcontent = []
                for line in content.splitlines():
                    # Match the pattern
                    matches_pattern = re.findall(rf"^.*{pattern[0]}.*$", line)
                    if matches_pattern:
                        found_match = True

                        # Replace the match
                        newline = re.sub(pattern[0], pattern[1], line)
                        newcontent.append(newline)
                    # No match, keep line as it is
                    else:
                        newcontent.append(line)

                if found_match:
                    content = "\n".join(newcontent) + "\n"
                else:
                    self.failed.append(
                        (f"Did not find pattern {pattern[0]} in module {module.component_name}", module.component_name)
                    )
                    return False

            # Write new content to the file
            with open(module.main_nf, "w") as fh:
                fh.write(content)

            # change version in environment.yml
            if not module.environment_yml:
                log.error(f"Could not read `environment.yml` of {module.component_name} module.")
                return False
            with open(module.environment_yml) as fh:
                env_yml = yaml.safe_load(fh)
            env_yml["dependencies"][0] = re.sub(
                bioconda_packages[0], f"bioconda::{bioconda_tool_name}={last_ver}", env_yml["dependencies"][0]
            )
            with open(module.environment_yml, "w") as fh:
                yaml.dump(env_yml, fh, default_flow_style=False, Dumper=custom_yaml_dumper())

            self.updated.append(
                (
                    f"Module updated:  {bioconda_version} --> {last_ver}",
                    module.component_name,
                )
            )
            return True

        else:
            self.up_to_date.append((f"Module version up to date: {module.component_name}", module.component_name))
            return True

    def get_bioconda_version(self, module: NFCoreComponent) -> list[str]:
        """
        Extract the bioconda version from a module
        """
        # Check whether file exists and load it
        bioconda_packages = []
        if module.environment_yml is not None and module.environment_yml.exists():
            with open(module.environment_yml) as fh:
                env_yml = yaml.safe_load(fh)
            bioconda_packages = env_yml.get("dependencies", [])
        else:
            log.error(f"Could not read `environment.yml` of {module.component_name} module.")

        return bioconda_packages

    def _print_results(self) -> None:
        """
        Print the results for the bump_versions command
        Uses the ``rich`` library to print a set of formatted tables to the command line
        summarising the linting results.
        """

        log.debug("Printing bump_versions results")

        console = Console(force_terminal=rich_force_colors())
        # Find maximum module name length
        max_mod_name_len = 40
        for m in [self.up_to_date, self.updated, self.failed]:
            try:
                max_mod_name_len = max(len(m[2]), max_mod_name_len)
            except Exception:
                pass

        def format_result(module_updates: list[tuple[str, str]], table: Table) -> Table:
            """
            Create rows for module updates
            """
            # TODO: Row styles don't work current as table-level style overrides.
            # I'd like to make an issue about this on the rich repo so leaving here in case there is a future fix
            last_modname = ""
            row_style = None
            for module_update in module_updates:
                if last_modname and module_update[1] != last_modname:
                    if row_style:
                        row_style = None
                    else:
                        row_style = "magenta"
                last_modname = module_update[1]
                table.add_row(
                    Markdown(f"{module_update[1]}"),
                    Markdown(f"{module_update[0]}"),
                    style=row_style,
                )
            return table

        # Table of up to date modules
        if len(self.up_to_date) > 0 and self.show_up_to_date:
            console.print(
                Panel(
                    rf"[!] {len(self.up_to_date)} Module{_s(self.up_to_date)} version{_s(self.up_to_date)} up to date.",
                    style="bold green",
                )
            )
            table = Table(style="green", box=ROUNDED)
            table.add_column("Module name", width=max_mod_name_len)
            table.add_column("Update Message")
            table = format_result(self.up_to_date, table)
            console.print(table)

        # Table of updated modules
        if len(self.updated) > 0:
            console.print(Panel(rf"[!] {len(self.updated)} Module{_s(self.updated)} updated", style="bold yellow"))
            table = Table(style="yellow", box=ROUNDED)
            table.add_column("Module name", width=max_mod_name_len)
            table.add_column("Update message")
            table = format_result(self.updated, table)
            console.print(table)

        # Table of modules that couldn't be updated
        if len(self.failed) > 0:
            console.print(Panel(rf"[!] {len(self.failed)} Module update{_s(self.failed)} failed", style="bold red"))
            table = Table(style="red", box=ROUNDED)
            table.add_column("Module name", width=max_mod_name_len)
            table.add_column("Update message")
            table = format_result(self.failed, table)
            console.print(table)

        # Table of modules ignored due to `.nf-core.yml`
        if len(self.ignored) > 0:
            console.print(Panel(rf"[!] {len(self.ignored)} Module update{_s(self.ignored)} ignored", style="grey58"))
            table = Table(style="grey58", box=ROUNDED)
            table.add_column("Module name", width=max_mod_name_len)
            table.add_column("Update message")
            table = format_result(self.ignored, table)
            console.print(table)
