import logging
import os
from pathlib import Path

import questionary
import yaml
from rich import box
from rich.console import Group
from rich.markdown import Markdown
from rich.panel import Panel
from rich.syntax import Syntax
from rich.table import Table
from rich.text import Text

import nf_core.utils
from nf_core.components.components_command import ComponentCommand
from nf_core.components.constants import NF_CORE_MODULES_REMOTE
from nf_core.modules.modules_json import ModulesJson

log = logging.getLogger(__name__)


class ComponentInfo(ComponentCommand):
    """
    Class to print information of a module/subworkflow.

    Attributes
    ----------
    meta : YAML object
        stores the information from meta.yml file
    local_path : str
        path of the local modules/subworkflows
    remote_location : str
        remote repository URL
    local : bool
        indicates if the module/subworkflow is locally installed or not
    repo_type : str
        repository type. Can be either 'pipeline' or 'modules'
    modules_json : ModulesJson object
        contains 'modules.json' file information from a pipeline
    component_name : str
        name of the module/subworkflow to get information from

    Methods
    -------
    init_mod_name(component)
        Makes sure that we have a modules/subworkflows name
    get_component_info()
        Given the name of a module/subworkflow, parse meta.yml and print usage help
    get_local_yaml()
        Attempt to get the meta.yml file from a locally installed module/subworkflow
    get_remote_yaml()
        Attempt to get the meta.yml file from a remote repo
    generate_component_info_help()
        Take the parsed meta.yml and generate rich help
    """

    def __init__(
        self,
        component_type: str,
        pipeline_dir: str | Path,
        component_name: str,
        remote_url: str | None = None,
        branch: str | None = None,
        no_pull: bool = False,
    ):
        super().__init__(component_type, pipeline_dir, remote_url, branch, no_pull)
        self.meta: dict | None = None
        self.local_path: Path | None = None
        self.remote_location: str | None = None
        self.local: bool = False
        self.modules_json: ModulesJson | None = None

        if self.repo_type == "pipeline":
            # Check modules directory structure
            if self.component_type == "modules":
                self.check_modules_structure()
            # Check modules.json up to date
            self.modules_json = ModulesJson(self.directory)
            self.modules_json.check_up_to_date()
        else:
            self.modules_json = None
        self.component = self.init_mod_name(component_name)

    def _configure_repo_and_paths(self, nf_dir_req=False) -> None:
        """
        Override the default with nf_dir_req set to False to allow
        info to be run from anywhere and still return remote info
        """
        return super()._configure_repo_and_paths(nf_dir_req)

    def init_mod_name(self, component: str | None) -> str:
        """
        Makes sure that we have a module/subworkflow name before proceeding.

        Args:
            module: str: Module name to check
        """
        if component is None:
            self.local = questionary.confirm(
                f"Is the {self.component_type[:-1]} locally installed?", style=nf_core.utils.nfcore_question_style
            ).unsafe_ask()
            if self.local:
                if self.repo_type == "modules":
                    components = self.get_components_clone_modules()
                elif self.repo_type == "pipeline":
                    assert self.modules_json is not None  # mypy
                    all_components: list[tuple[str, str]] = self.modules_json.get_all_components(
                        self.component_type
                    ).get(self.modules_repo.remote_url, [])

                    components = [
                        component if directory == self.modules_repo.repo_path else f"{directory}/{component}"
                        for directory, component in all_components
                    ]
                    if not components:
                        raise UserWarning(
                            f"No {self.component_type[:-1]} installed from '{self.modules_repo.remote_url}'"
                        )
                else:
                    raise UserWarning("Unknown repository type")
            else:
                components = self.modules_repo.get_avail_components(self.component_type)
            components.sort()
            component = questionary.autocomplete(
                f"Please select a {self.component_type[:-1]}",
                choices=components,
                style=nf_core.utils.nfcore_question_style,
            ).unsafe_ask()
            while component not in components:
                log.info(f"'{component}' is not a valid {self.component_type[:-1]} name")
                component = questionary.autocomplete(
                    f"Please select a new {self.component_type[:-1]}",
                    choices=components,
                    style=nf_core.utils.nfcore_question_style,
                ).unsafe_ask()
        else:
            if self.repo_type == "pipeline":
                assert self.modules_json is not None  # mypy
                # check if the module is locally installed
                local_paths = self.modules_json.get_all_components(self.component_type).get(
                    self.modules_repo.remote_url
                )  # type: ignore
                if local_paths is not None:
                    for directory, comp in local_paths:
                        if comp == component:
                            component_base_path = Path(self.directory, self.component_type)
                            self.local_path = Path(component_base_path, directory, component)
                            break
                if self.local_path:
                    self.local = True

        return component

    def get_component_info(self):
        """Given the name of a module/subworkflow, parse meta.yml and print usage help."""

        # Running with a local install, try to find the local meta
        if self.local:
            self.meta = self.get_local_yaml()

        # Either failed locally or in remote mode
        if not self.meta:
            self.meta = self.get_remote_yaml()

        # Could not find the meta
        if self.meta is None:
            raise UserWarning(f"Could not find {self.component_type[:-1]} '{self.component}'")

        return self.generate_component_info_help()

    def get_local_yaml(self) -> dict | None:
        """Attempt to get the meta.yml file from a locally installed module/subworkflow.

        Returns:
            dict | None: Parsed meta.yml if found, None otherwise
        """

        if self.repo_type == "pipeline":
            assert self.modules_json is not None  # mypy
            # Try to find and load the meta.yml file
            component_base_path = Path(self.directory, self.component_type)
            # Check that we have any modules/subworkflows installed from this repo
            components = self.modules_json.get_all_components(self.component_type).get(self.modules_repo.remote_url)
            if components is None:
                raise LookupError(f"No {self.component_type[:-1]} installed from {self.modules_repo.remote_url}")

            component_names = [component for _, component in components]

            if self.component in component_names:
                install_dir = [directory for directory, module in components if module == self.component][0]
                comp_dir = Path(component_base_path, install_dir, self.component)
                meta_fn = Path(comp_dir, "meta.yml")
                if meta_fn.exists():
                    log.debug(f"Found local file: {meta_fn}")
                    with open(meta_fn) as fh:
                        self.local_path = comp_dir
                        return yaml.safe_load(fh)

            log.debug(f"{self.component_type[:-1].title()} '{self.component}' meta.yml not found locally")
        else:
            component_base_path = Path(self.directory, self.component_type, self.org)
            if self.component in os.listdir(component_base_path):
                comp_dir = Path(component_base_path, self.component)
                meta_fn = Path(comp_dir, "meta.yml")
                if meta_fn.exists():
                    log.debug(f"Found local file: {meta_fn}")
                    with open(meta_fn) as fh:
                        self.local_path = comp_dir
                        return yaml.safe_load(fh)
            log.debug(f"{self.component_type[:-1].title()} '{self.component}' meta.yml not found locally")

        return {}

    def get_remote_yaml(self) -> dict | None:
        """Attempt to get the meta.yml file from a remote repo.

        Returns:
            dict or bool: Parsed meta.yml found, False otherwise
        """
        # Check if our requested module/subworkflow is there
        if self.component not in self.modules_repo.get_avail_components(self.component_type):
            return None

        file_contents = self.modules_repo.get_meta_yml(self.component_type, self.component)
        if file_contents is None:
            return None
        self.remote_location = self.modules_repo.remote_url
        return yaml.safe_load(file_contents)

    def generate_params_table(self, type) -> Table:
        "Generate a rich table for inputs and outputs"
        table = Table(expand=True, show_lines=True, box=box.MINIMAL_HEAVY_HEAD, padding=0)
        table.add_column(f":inbox_tray: {type}")
        table.add_column("Description")
        if self.component_type == "modules":
            table.add_column("Pattern", justify="right", style="green")
        elif self.component_type == "subworkflows":
            table.add_column("Structure", justify="right", style="green")
        return table

    def get_channel_structure(self, structure: dict) -> str:
        "Get the structure of a channel"
        structure_str = ""
        for key, info in structure.items():
            pattern = f" - {info['pattern']}" if info.get("pattern") else ""
            structure_str += f"{key} ({info['type']}{pattern})"
        return structure_str

    def generate_component_info_help(self):
        """Take the parsed meta.yml and generate rich help.

        Returns:
            rich renderable
        """

        renderables = []

        # Intro panel
        intro_text = Text()
        if self.local_path:
            intro_text.append(Text.from_markup(f"Location: [blue]{self.local_path}\n"))
        elif self.remote_location:
            intro_text.append(
                Text.from_markup(
                    ":globe_with_meridians: Repository: "
                    f"{'[link={self.remote_location}]' if self.remote_location.startswith('http') else ''}"
                    f"{self.remote_location}"
                    f"{'[/link]' if self.remote_location.startswith('http') else ''}"
                    "\n"
                )
            )
        if self.meta is None:
            raise UserWarning("No meta.yml file found")
        if self.meta.get("tools"):
            tools_strings = []
            for tool in self.meta["tools"]:
                for tool_name, tool_meta in tool.items():
                    if "homepage" in tool_meta:
                        tools_strings.append(f"[link={tool_meta['homepage']}]{tool_name}[/link]")
                    else:
                        tools_strings.append(f"{tool_name}")
            intro_text.append(Text.from_markup(f":wrench: Tools: {', '.join(tools_strings)}\n", style="dim"))

        if self.meta.get("description"):
            intro_text.append(Text.from_markup(f":book: Description: {self.meta['description']}", style="dim"))

        renderables.append(
            Panel(
                intro_text,
                title=f"[bold]{self.component_type[:-1].title()}: [green]{self.component}\n",
                title_align="left",
            )
        )

        # Inputs
        if self.meta.get("input"):
            inputs_table = self.generate_params_table("Inputs")
            for i, input in enumerate(self.meta["input"]):
                inputs_table.add_row(f"[italic]input[{i}][/]", "", "")
                if self.component_type == "modules":
                    for element in input:
                        for key, info in element.items():
                            inputs_table.add_row(
                                f"[orange1 on black] {key} [/][dim i] ({info['type']})",
                                Markdown(info["description"] if info["description"] else ""),
                                info.get("pattern", ""),
                            )
                elif self.component_type == "subworkflows":
                    for key, info in input.items():
                        inputs_table.add_row(
                            f"[orange1 on black] {key} [/][dim i]",
                            Markdown(info["description"] if info["description"] else ""),
                            self.get_channel_structure(info["structure"]) if info.get("structure") else "",
                        )

            renderables.append(inputs_table)

        # Outputs
        if self.meta.get("output"):
            outputs_table = self.generate_params_table("Outputs")
            if self.component_type == "modules":
                for ch_name, elements in self.meta["output"].items():
                    outputs_table.add_row(f"{ch_name}", "", "")
                    for element in elements:
                        if isinstance(element, list):
                            for ch in element:
                                for key, info in ch.items():
                                    outputs_table.add_row(
                                        f"[orange1 on black] {key} [/][dim i] ({info['type']})",
                                        Markdown(info["description"] if info["description"] else ""),
                                        info.get("pattern", ""),
                                    )
                        elif isinstance(element, dict):
                            for key, info in element.items():
                                outputs_table.add_row(
                                    f"[orange1 on black] {key} [/][dim i] ({info['type']})",
                                    Markdown(info["description"] if info["description"] else ""),
                                    info.get("pattern", ""),
                                )
            elif self.component_type == "subworkflows":
                for output in self.meta["output"]:
                    for key, info in output.items():
                        outputs_table.add_row(
                            f"[orange1 on black] {key} [/][dim i]",
                            Markdown(info["description"] if info["description"] else ""),
                            self.get_channel_structure(info["structure"]) if info.get("structure") else "",
                        )

            renderables.append(outputs_table)

        # Installation command
        if self.remote_location and not self.local:
            cmd_base = f"nf-core {self.component_type}"
            if self.remote_location != NF_CORE_MODULES_REMOTE:
                cmd_base = f"nf-core {self.component_type} --git-remote {self.remote_location}"
            renderables.append(
                Text.from_markup(f"\n :computer:  Installation command: [magenta]{cmd_base} install {self.component}\n")
            )

        # Print include statement
        if self.local_path and self.modules_repo.repo_path is not None:
            install_folder = Path(self.directory, self.component_type, self.modules_repo.repo_path)
            component_name = "_".join(self.component.upper().split("/"))
            renderables.append(
                Text.from_markup(f"\n [blue]Use the following statement to include this {self.component_type[:-1]}:")
            )
            renderables.append(
                Syntax(
                    f"include {{ {component_name} }} from '../{Path(install_folder, self.component).relative_to(self.directory)}/main'",
                    "groovy",
                    theme="ansi_dark",
                    padding=1,
                )
            )
            if self.component_type == "subworkflows":
                subworkflow_config = Path(install_folder, self.component, "nextflow.config").relative_to(self.directory)
                if os.path.isfile(subworkflow_config):
                    renderables.append(
                        Text.from_markup("\n [blue]Add the following config statement to use this subworkflow:")
                    )
                    renderables.append(
                        Syntax(f"includeConfig '{subworkflow_config}'", "groovy", theme="ansi_dark", padding=1)
                    )

        return Group(*renderables)
