import logging
from pathlib import Path

import questionary
from rich.console import Console, Group
from rich.panel import Panel
from rich.syntax import Syntax

import nf_core.utils
from nf_core.components.components_command import ComponentCommand
from nf_core.modules.modules_json import ModulesJson

from .install import ComponentInstall

log = logging.getLogger(__name__)


class ComponentRemove(ComponentCommand):
    def __init__(self, component_type, pipeline_dir, remote_url=None, branch=None, no_pull=False):
        super().__init__(component_type, pipeline_dir, remote_url, branch, no_pull)

    def remove(self, component, repo_url=None, repo_path=None, removed_by=None, removed_components=None, force=False):
        """
        Remove an already installed module/subworkflow
        This command only works for modules/subworkflows that are installed from 'nf-core/modules'

        Args:
            component (str): Name of the component to remove
            repo_url (str): URL of the repository where the component is located
            repo_path (str): Directory where the component is installed
            removed_by (str): Name of the component that is removing the current component
                (a subworkflow name if the component is a dependency or "modules" or "subworkflows" if it is not a dependency)
            removed_components (list[str]): list of components that have been removed during a recursive remove of subworkflows
            force (bool): Force removal of component, even if there is still an include statement in a workflow file

        Returns:
            bool: True if any item has been removed, False if not
        """
        if self.repo_type == "modules":
            log.error(f"You cannot remove a {self.component_type[:-1]} in a clone of nf-core/modules")
            return False

        # Check modules directory structure
        if self.component_type == "modules":
            self.check_modules_structure()

        # Check whether pipeline is valid and with a modules.json file
        self.has_valid_directory()
        self.has_modules_file()

        if repo_path is None:
            repo_path = self.modules_repo.repo_path
        if repo_url is None:
            repo_url = self.modules_repo.remote_url
        if component is None:
            component = questionary.autocomplete(
                f"{self.component_type[:-1]} name:",
                choices=self.components_from_repo(repo_path),
                style=nf_core.utils.nfcore_question_style,
            ).unsafe_ask()

        if removed_components is None:
            removed_components = []

        # Get the module/subworkflow directory
        component_dir = Path(self.directory, self.component_type, repo_path, component)

        # Load the modules.json file
        modules_json = ModulesJson(self.directory)
        modules_json.load()

        # Verify that the module/subworkflow is actually installed
        if not component_dir.exists():
            log.error(f"Installation directory '{component_dir}' does not exist.")

            if modules_json.component_present(component, repo_url, repo_path, self.component_type):
                log.error(f"Found entry for '{component}' in 'modules.json'. Removing...")
                modules_json.remove_entry(self.component_type, component, repo_url, repo_path)
            return False

        # remove all dependent components based on installed_by entry
        # Remove entry from modules.json
        removed = False
        removed_components = []
        # Remove component from modules.json
        removed_component = modules_json.remove_entry(
            self.component_type,
            component,
            repo_url,
            repo_path,
            removed_by=removed_by,
        )
        removed_component_dir = Path(self.component_type, repo_path, component)
        if removed_component:
            # check if the module/subworkflow has been manually included in the pipeline
            include_stmts = self.check_if_in_include_stmts(str(removed_component_dir))
            if include_stmts:
                # print the include statements
                log.warning(
                    f"The {self.component_type[:-1]} '{component}' is still included in the following workflow file{nf_core.utils.plural_s(include_stmts)}:"
                )
                console = Console()
                for file, stmts in include_stmts.items():
                    renderables = []
                    for stmt in stmts:
                        # check that the line number is integer
                        if not isinstance(stmt["line_number"], int):
                            log.error(
                                f"Could not parse line number '{stmt['line_number']}' in '{file}'. Please report this issue."
                            )
                            continue

                        renderables.append(
                            Syntax(
                                str(stmt["line"]),
                                "groovy",
                                theme="ansi_dark",
                                line_numbers=True,
                                start_line=stmt["line_number"],
                            )
                        )
                    console.print(
                        Panel(
                            Group(*renderables),
                            title=f"{file}",
                            style="white",
                            title_align="center",
                            padding=1,
                        )
                    )
                # ask the user if they still want to remove the component, add it back otherwise
                if not force:
                    if not questionary.confirm(
                        f"Do you still want to remove the {self.component_type[:-1]} '{component}'?",
                        style=nf_core.utils.nfcore_question_style,
                    ).unsafe_ask():
                        # add the component back to modules.json
                        if not ComponentInstall(self.directory, self.component_type, force=True).install(
                            component, silent=True
                        ):
                            log.warning(
                                f"Could not install the {self.component_type[:-1]} '{component}', please install it manually with 'nf-core {self.component_type} install  {component}'."
                            )
                        removed_components.append(component)
                        return removed
            # Remove the component files of all entries removed from modules.json
            removed = (
                True
                if self.clear_component_dir(component, Path(self.directory, removed_component_dir)) or removed
                else False
            )
            removed_components.append(component)

        if removed:
            if self.component_type == "subworkflows":
                removed_by = component
                dependent_components = modules_json.get_dependent_components(self.component_type, component, {})
                for component_name, component_data in dependent_components.items():
                    component_repo, component_install_dir, component_type = component_data
                    if component_name in removed_components:
                        continue
                    original_component_type = self.component_type
                    self.component_type = component_type
                    dependency_removed = self.remove(
                        component_name,
                        component_repo,
                        component_install_dir,
                        removed_by=removed_by,
                        removed_components=removed_components,
                    )
                    self.component_type = original_component_type
                    # remember removed dependencies
                    if dependency_removed:
                        removed_components.append(component_name.replace("/", "_"))
            # print removed dependencies
            if removed_components:
                log.info(f"Removed files for '{component}' and its dependencies '{', '.join(removed_components)}'.")
            else:
                log.info(f"Removed files for '{component}'.")
        else:
            installed_by = modules_json.get_installed_by_entries(self.component_type, component)
            if installed_by == [self.component_type]:
                log.error(
                    f"Did not remove '{component}', because it was also manually installed. Only updated 'installed_by' entry in modules.json."
                )
            log.info(
                f"""Did not remove {self.component_type[:-1]} '{component}', because it was also installed by {", ".join(f"'{d}'" for d in installed_by)}. Only updated the 'installed_by' entry in modules.json."""
            )
        return removed
