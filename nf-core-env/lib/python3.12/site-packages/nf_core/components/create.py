"""
The ComponentCreate class handles generating of module and subworkflow templates
"""

import glob
import json
import logging
import re
import shutil
import subprocess
from pathlib import Path

import jinja2
import questionary
import rich
import rich.prompt
import ruamel.yaml
from packaging.version import parse as parse_version

import nf_core
import nf_core.utils
from nf_core.components.components_command import ComponentCommand
from nf_core.components.components_utils import get_biotools_id, get_biotools_response, get_channel_info_from_biotools
from nf_core.pipelines.lint_utils import run_prettier_on_file

log = logging.getLogger(__name__)

# Set yaml options for meta.yml files
ruamel.yaml.representer.RoundTripRepresenter.ignore_aliases = (
    lambda x, y: True
)  # Fix to not print aliases. https://stackoverflow.com/a/64717341
yaml = ruamel.yaml.YAML()
yaml.preserve_quotes = True
yaml.indent(mapping=2, sequence=2, offset=0)


class ComponentCreate(ComponentCommand):
    def __init__(
        self,
        component_type: str,
        directory: Path = Path("."),
        component: str = "",
        author: str | None = None,
        process_label: str | None = None,
        has_meta: str | None = None,
        force: bool = False,
        conda_name: str | None = None,
        conda_version: str | None = None,
        empty_template: bool = False,
        migrate_pytest: bool = False,  # TODO: Deprecate this flag in the future
    ):
        super().__init__(component_type, directory)
        self.directory = directory
        self.component = component
        self.author = author
        self.process_label = process_label
        self.has_meta = has_meta
        self.force_overwrite = force
        self.subtool = None
        self.tool_conda_name = conda_name
        self.tool_conda_version = conda_version
        self.tool_licence = ""
        self.tool_description = ""
        self.tool_doc_url = ""
        self.tool_dev_url = ""
        self.bioconda = None
        self.singularity_container = None
        self.docker_container = None
        self.file_paths: dict[str, Path] = {}
        self.not_empty_template = not empty_template
        self.migrate_pytest = migrate_pytest
        self.tool_identifier = ""

    def create(self) -> bool:
        """
        Create a new DSL2 module or subworkflow from the nf-core template.

        A module should be named just <tool> or <tool/subtool>
        e.g fastqc or samtools/sort, respectively.

        The subworkflow should be named as the main file type it operates on and a short description of the task performed
        e.g bam_sort or bam_sort_samtools, respectively.

        If <directory> is a pipeline, this function creates a file called:
        '<directory>/modules/local/tool/main.nf'
            OR
        '<directory>/modules/local/tool/subtool/main.nf'
            OR for subworkflows
        '<directory>/subworkflows/local/subworkflow_name/main.nf'

        If <directory> is a clone of nf-core/modules, it creates or modifies the following files:

        For modules:

        ```tree
        modules/nf-core/tool/subtool/
        ├── main.nf
        ├── meta.yml
        ├── environment.yml
        └── tests
            └── main.nf.test
        ```

        The function will attempt to automatically find a Bioconda package called <component>
        and matching Docker / Singularity images from BioContainers.

        For subworkflows:

        ```tree
        subworkflows/nf-core/tool/subtool/
        ├── main.nf
        ├── meta.yml
        └── tests
            └── main.nf.test
        ```

        """
        if self.component_type == "modules":
            # Check modules directory structure
            self.check_modules_structure()

        # Check whether the given directory is a nf-core pipeline or a clone of nf-core/modules

        log.info(f"Repository type: [blue]{self.repo_type}")
        if self.directory != ".":
            log.info(f"Base directory: '{self.directory}'")

        log.info(
            "[yellow]Press enter to use default values [cyan bold](shown in brackets)[/] [yellow]or type your own responses. "
            "ctrl+click [link=https://youtu.be/dQw4w9WgXcQ]underlined text[/link] to open links."
        )

        # Collect component info via prompt if empty or invalid
        self._collect_name_prompt()

        # Determine the component name
        self.component_name = self.component
        self.component_dir = Path(self.component)

        if self.subtool:
            self.component_name = f"{self.component}/{self.subtool}"
            self.component_dir = Path(self.component, self.subtool)

        self.component_name_underscore = self.component_name.replace("/", "_")

        # Check existence of directories early for fast-fail
        self.file_paths = self._get_component_dirs()

        if self.migrate_pytest:
            # Rename the component directory to old
            component_old_dir = Path(str(self.component_dir) + "_old")
            component_parent_path = Path(self.directory, self.component_type, self.org)
            component_old_path = component_parent_path / component_old_dir
            component_path = component_parent_path / self.component_dir

            component_path.rename(component_old_path)
        else:
            if self.component_type == "modules":
                # Try to find a bioconda package for 'component'
                self._get_bioconda_tool()
                name = self.tool_conda_name if self.tool_conda_name else self.component
                # Try to find a biotools entry for 'component'
                biotools_data = get_biotools_response(name)
                if biotools_data:
                    self.tool_identifier = get_biotools_id(biotools_data, name)
                    # Obtain EDAM ontologies for inputs and outputs
                    channel_info = get_channel_info_from_biotools(biotools_data, name)
                    if channel_info:
                        self.inputs, self.outputs = channel_info

            # Prompt for GitHub username
            self._get_username()

            if self.component_type == "modules":
                self._get_module_structure_components()

        # Add a valid organization name for nf-test tags
        not_alphabet = re.compile(r"[^a-zA-Z]")
        self.org_alphabet = not_alphabet.sub("", self.org)

        # Create component template with jinja2
        assert self._render_template()
        log.info(f"Created component template: '{self.component_name}'")

        if self.component_type == "modules":
            # Generate meta.yml inputs and outputs
            self.generate_meta_yml_file()

        if self.migrate_pytest:
            self._copy_old_files(component_old_path)
            log.info("Migrate pytest tests: Copied original module files to new module")
            shutil.rmtree(component_old_path)
            self._print_and_delete_pytest_files()

        new_files = [str(path) for path in self.file_paths.values()]

        run_prettier_on_file(new_files)

        log.info("Created following files:\n  " + "\n  ".join(new_files))
        return True

    def _get_bioconda_tool(self):
        """
        Try to find a bioconda package for 'tool'
        """
        while True:
            try:
                if self.tool_conda_name:
                    anaconda_response = nf_core.utils.anaconda_package(self.tool_conda_name, ["bioconda"])
                else:
                    anaconda_response = nf_core.utils.anaconda_package(self.component, ["bioconda"])

                if not self.tool_conda_version:
                    version = anaconda_response.get("latest_version")
                    if not version:
                        version = str(max([parse_version(v) for v in anaconda_response["versions"]]))
                else:
                    version = self.tool_conda_version

                self.tool_licence = nf_core.utils.parse_anaconda_licence(anaconda_response, version)
                self.tool_description = anaconda_response.get("summary", "")
                self.tool_doc_url = anaconda_response.get("doc_url", "")
                self.tool_dev_url = anaconda_response.get("dev_url", "")
                if self.tool_conda_name:
                    self.bioconda = "bioconda::" + self.tool_conda_name + "=" + version
                else:
                    self.bioconda = "bioconda::" + self.component + "=" + version
                log.info(f"Using Bioconda package: '{self.bioconda}'")
                break
            except (ValueError, LookupError) as e:
                log.warning(
                    f"Could not find Conda dependency using the Anaconda API: '{self.tool_conda_name if self.tool_conda_name else self.component}'"
                )
                if rich.prompt.Confirm.ask("[violet]Do you want to enter a different Bioconda package name?"):
                    self.tool_conda_name = rich.prompt.Prompt.ask("[violet]Name of Bioconda package").strip()
                    continue
                else:
                    log.warning(
                        f"{e}\nBuilding module without tool software and meta, you will need to enter this information manually."
                    )
                    break

        # Try to get the container tag (only if bioconda package was found)
        if self.bioconda:
            try:
                if self.tool_conda_name:
                    self.docker_container, self.singularity_container = nf_core.utils.get_biocontainer_tag(
                        self.tool_conda_name, version
                    )
                else:
                    self.docker_container, self.singularity_container = nf_core.utils.get_biocontainer_tag(
                        self.component, version
                    )
                log.info(f"Using Docker container: '{self.docker_container}'")
                log.info(f"Using Singularity container: '{self.singularity_container}'")
            except (ValueError, LookupError) as e:
                log.info(f"Could not find a Docker/Singularity container ({e})")

    def _get_module_structure_components(self):
        process_label_defaults = [
            "process_single",
            "process_low",
            "process_medium",
            "process_high",
            "process_long",
            "process_high_memory",
        ]
        if self.process_label is None:
            log.info(
                "Provide an appropriate resource label for the process, taken from the "
                "[link=https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/conf/base.config#L29]nf-core pipeline template[/link].\n"
                "For example: {}".format(", ".join(process_label_defaults))
            )
        while self.process_label is None:
            self.process_label = questionary.autocomplete(
                "Process resource label:",
                choices=process_label_defaults,
                style=nf_core.utils.nfcore_question_style,
                default="process_single",
            ).unsafe_ask()

        if self.has_meta is None:
            log.info(
                "Where applicable all sample-specific information e.g. 'id', 'single_end', 'read_group' "
                "MUST be provided as an input via a Groovy Map called 'meta'. "
                "This information may [italic]not[/] be required in some instances, for example "
                "[link=https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf]indexing reference genome files[/link]."
            )
        while self.has_meta is None:
            self.has_meta = rich.prompt.Confirm.ask(
                "[violet]Will the module require a meta map of sample information?",
                default=True,
            )

    def _render_template(self) -> bool | None:
        """
        Create new module/subworkflow files with Jinja2.
        """
        object_attrs = vars(self)
        # Run jinja2 for each file in the template folder
        env = jinja2.Environment(
            loader=jinja2.PackageLoader("nf_core", f"{self.component_type[:-1]}-template"),
            keep_trailing_newline=True,
        )
        for template_fn, dest_fn in self.file_paths.items():
            log.debug(f"Rendering template file: '{template_fn}'")
            j_template = env.get_template(template_fn)
            object_attrs["nf_core_version"] = nf_core.__version__
            try:
                rendered_output = j_template.render(object_attrs)
            except Exception as e:
                log.error(f"Could not render template file '{template_fn}':\n{e}")
                raise e

            # Write output to the target file
            log.debug(f"Writing output to: '{dest_fn}'")
            dest_fn.parent.mkdir(exist_ok=True, parents=True)
            with open(dest_fn, "w") as fh:
                log.debug(f"Writing output to: '{dest_fn}'")
                fh.write(rendered_output)

            # Mirror file permissions
            template_stat = (
                Path(nf_core.__file__).parent / f"{self.component_type[:-1]}-template" / template_fn
            ).stat()
            dest_fn.chmod(template_stat.st_mode)
        return True

    def _collect_name_prompt(self):
        """
        Collect module/subworkflow info via prompt if empty or invalid
        """
        # Collect module info via prompt if empty or invalid
        self.subtool = None
        if self.component_type == "modules":
            pattern = r"[^a-z\d/]"
        elif self.component_type == "subworkflows":
            pattern = r"[^a-z\d_/]"
        if self.component is None:
            self.component = ""
        while self.component == "" or re.search(pattern, self.component) or self.component.count("/") > 0:
            # Check + auto-fix for invalid chacters
            if re.search(pattern, self.component):
                if self.component_type == "modules":
                    log.warning("Tool/subtool name must be lower-case letters only, with no punctuation")
                elif self.component_type == "subworkflows":
                    log.warning("Subworkflow name must be lower-case letters only, with no punctuation")
                name_clean = re.sub(r"[^a-z\d/]", "", self.component.lower())
                if rich.prompt.Confirm.ask(f"[violet]Change '{self.component}' to '{name_clean}'?"):
                    self.component = name_clean
                else:
                    self.component = ""

            if self.component_type == "modules":
                # Split into tool and subtool
                if self.component.count("/") > 1:
                    log.warning("Tool/subtool can have maximum one '/' character")
                    self.component = ""
                elif self.component.count("/") == 1:
                    self.component, self.subtool = self.component.split("/")
                else:
                    self.subtool = None  # Reset edge case: entered '/subtool' as name and gone round loop again

            # Prompt for new entry if we reset
            if self.component == "":
                if self.component_type == "modules":
                    self.component = rich.prompt.Prompt.ask("[violet]Name of tool/subtool").strip()
                elif self.component_type == "subworkflows":
                    self.component = rich.prompt.Prompt.ask("[violet]Name of subworkflow").strip()

    def _get_component_dirs(self) -> dict[str, Path]:
        """Given a directory and a tool/subtool or subworkflow, set the file paths and check if they already exist

        Returns dict: keys are relative paths to template files, vals are target paths.
        """
        file_paths = {}
        if self.repo_type == "pipeline":
            component_dir = Path(self.directory, self.component_type, "local", self.component_dir)

        elif self.repo_type == "modules":
            component_dir = Path(self.directory, self.component_type, self.org, self.component_dir)
        else:
            raise ValueError("`repo_type` not set correctly")

        # Check if module/subworkflow directories exist already
        if component_dir.exists() and not self.force_overwrite and not self.migrate_pytest:
            raise UserWarning(
                f"{self.component_type[:-1]} directory exists: '{component_dir}'. Use '--force' to overwrite"
            )

        if self.component_type == "modules":
            # If a subtool, check if there is a module called the base tool name already
            parent_tool_main_nf = Path(
                self.directory,
                self.component_type,
                self.org,
                self.component,
                "main.nf",
            )
            if self.subtool and parent_tool_main_nf.exists() and not self.migrate_pytest:
                raise UserWarning(
                    f"Module '{parent_tool_main_nf}' exists already, cannot make subtool '{self.component_name}'"
                )

            # If no subtool, check that there isn't already a tool/subtool
            tool_glob = glob.glob(f"{Path(self.directory, self.component_type, self.org, self.component)}/*/main.nf")
            if not self.subtool and tool_glob and not self.migrate_pytest:
                raise UserWarning(
                    f"Module subtool '{tool_glob[0]}' exists already, cannot make tool '{self.component_name}'"
                )
        # Set file paths
        # For modules - can be tool/ or tool/subtool/ so can't do in template directory structure
        file_paths["main.nf"] = component_dir / "main.nf"
        file_paths["meta.yml"] = component_dir / "meta.yml"
        if self.component_type == "modules":
            file_paths["environment.yml"] = component_dir / "environment.yml"
        file_paths["tests/main.nf.test.j2"] = component_dir / "tests" / "main.nf.test"

        return file_paths

    def _get_username(self):
        """
        Prompt for GitHub username
        """
        # Try to guess the current user if `gh` is installed
        author_default = None
        try:
            gh_auth_user = json.loads(subprocess.check_output(["gh", "api", "/user"], stderr=subprocess.DEVNULL))
            author_default = f"@{gh_auth_user['login']}"
        except Exception as e:
            log.debug(f"Could not find GitHub username using 'gh' cli command: [red]{e}")

        # Regex to valid GitHub username: https://github.com/shinnn/github-username-regex
        github_username_regex = re.compile(r"^@[a-zA-Z\d](?:[a-zA-Z\d]|-(?=[a-zA-Z\d])){0,38}$")
        while self.author is None or not github_username_regex.match(self.author):
            if self.author is not None and not github_username_regex.match(self.author):
                log.warning("Does not look like a valid GitHub username (must start with an '@')!")
            self.author = rich.prompt.Prompt.ask(
                f"[violet]GitHub Username:[/]{' (@author)' if author_default is None else ''}",
                default=author_default,
            )

    def _copy_old_files(self, component_old_path):
        """Copy files from old module to new module"""
        log.debug("Copying original main.nf file")
        shutil.copyfile(component_old_path / "main.nf", self.file_paths["main.nf"])
        log.debug("Copying original meta.yml file")
        shutil.copyfile(component_old_path / "meta.yml", self.file_paths["meta.yml"])
        if self.component_type == "modules":
            log.debug("Copying original environment.yml file")
            shutil.copyfile(
                component_old_path / "environment.yml",
                self.file_paths["environment.yml"],
            )
            if (component_old_path / "templates").is_dir():
                log.debug("Copying original templates directory")
                shutil.copytree(
                    component_old_path / "templates",
                    self.file_paths["environment.yml"].parent / "templates",
                )
        # Create a nextflow.config file if it contains information other than publishDir
        pytest_dir = Path(self.directory, "tests", self.component_type, self.org, self.component_dir)
        nextflow_config = pytest_dir / "nextflow.config"
        if nextflow_config.is_file():
            with open(nextflow_config) as fh:
                config_lines = ""
                for line in fh:
                    if "publishDir" not in line and line.strip() != "":
                        config_lines += line
            # if the nextflow.config file only contained publishDir, non_publish_dir_lines will be 11 characters long (`process {\n}`)
            if len(config_lines) > 11:
                log.debug("Copying nextflow.config file from pytest tests")
                with open(
                    Path(
                        self.directory,
                        self.component_type,
                        self.org,
                        self.component_dir,
                        "tests",
                        "nextflow.config",
                    ),
                    "w+",
                ) as ofh:
                    ofh.write(config_lines)

    def _print_and_delete_pytest_files(self):
        """Prompt if pytest files should be deleted and printed to stdout"""
        pytest_dir = Path(self.directory, "tests", self.component_type, self.org, self.component_dir)
        if rich.prompt.Confirm.ask(
            "[violet]Do you want to delete the pytest files?[/]\nPytest file 'main.nf' will be printed to standard output to allow migrating the tests manually to 'main.nf.test'.",
            default=False,
        ):
            with open(pytest_dir / "main.nf") as fh:
                log.info(fh.read())
            if pytest_dir.is_symlink():
                resolved_dir = pytest_dir.resolve()
                log.debug(f"Removing symlink: {resolved_dir}")
                shutil.rmtree(resolved_dir)
                pytest_dir.unlink()
            else:
                shutil.rmtree(pytest_dir)
            log.info(
                "[yellow]Please convert the pytest tests to nf-test in 'main.nf.test'.[/]\n"
                "You can find more information about nf-test [link=https://nf-co.re/docs/contributing/modules#migrating-from-pytest-to-nf-test]at the nf-core web[/link]. "
            )
        else:
            log.info(
                "[yellow]Please migrate the pytest tests to nf-test in 'main.nf.test'.[/]\n"
                "You can find more information about nf-test [link=https://nf-co.re/docs/contributing/modules#migrating-from-pytest-to-nf-test]at the nf-core web[/link].\n"
                f"Once done, make sure to delete the module pytest files to avoid linting errors: {pytest_dir}"
            )
        # Delete tags from pytest_modules.yml
        modules_yml = Path(self.directory, "tests", "config", "pytest_modules.yml")
        with open(modules_yml) as fh:
            yml_file = yaml.load(fh)
        yml_key = str(self.component_dir) if self.component_type == "modules" else f"subworkflows/{self.component_dir}"
        if yml_key in yml_file:
            del yml_file[yml_key]
        with open(modules_yml, "w") as fh:
            yaml.dump(yml_file, fh)
        run_prettier_on_file(modules_yml)

    def generate_meta_yml_file(self) -> None:
        """
        Generate the meta.yml file.
        """
        # TODO: The meta.yml could be handled with a Pydantic model. The reason it is not implemented is because we want to maintain comments in the meta.yml file.
        with open(self.file_paths["meta.yml"]) as fh:
            meta_yml: ruamel.yaml.comments.CommentedMap = yaml.load(fh)

        versions: dict[str, list | dict] = {
            f"versions_{self.component}": [
                [
                    {"${task.process}": {"type": "string", "description": "The name of the process"}},
                    {f"{self.component}": {"type": "string", "description": "The name of the tool"}},
                    {
                        f"{self.component} --version": {
                            "type": "eval",
                            "description": "The expression to obtain the version of the tool",
                        },
                    },
                ]
            ]
        }

        versions_topic: dict[str, list | dict] = {
            "versions": [
                [
                    {"${task.process}": {"type": "string", "description": "The name of the process"}},
                    {f"{self.component}": {"type": "string", "description": "The name of the tool"}},
                    {
                        f"{self.component} --version": {
                            "type": "eval",
                            "description": "The expression to obtain the version of the tool",
                        },
                    },
                ]
            ]
        }

        if self.not_empty_template:
            meta_yml.yaml_set_comment_before_after_key(
                "name", before="# TODO nf-core: Add a description of the module and list keywords"
            )
            meta_yml["tools"][0].yaml_set_start_comment(
                "## TODO nf-core: Add a description and other details for the software below"
            )
            meta_yml["input"].yaml_set_start_comment(
                "### TODO nf-core: Add a description of all of the variables used as input", indent=2
            )
            meta_yml["output"].yaml_set_start_comment(
                "### TODO nf-core: Add a description of all of the variables used as output", indent=2
            )
            meta_yml["topics"].yaml_set_start_comment(
                "### TODO nf-core: Add a description of all of the variables used as topics", indent=2
            )

            if hasattr(self, "inputs") and len(self.inputs) > 0:
                inputs_array: list[dict | list[dict]] = []
                for i, (input_name, ontologies) in enumerate(self.inputs.items()):
                    channel_entry: dict[str, dict] = {
                        input_name: {
                            "type": "file",
                            "description": f"{input_name} file",
                            "pattern": f"*.{{{','.join(ontologies[2])}}}",
                            "ontologies": [
                                ruamel.yaml.comments.CommentedMap({"edam": f"{ont_id}"}) for ont_id in ontologies[0]
                            ],
                        }
                    }
                    for j, ont_desc in enumerate(ontologies[1]):
                        channel_entry[input_name]["ontologies"][j].yaml_add_eol_comment(ont_desc, "edam")
                    if self.has_meta:
                        meta_suffix = str(i + 1) if i > 0 else ""
                        meta_entry: dict[str, dict] = {
                            f"meta{meta_suffix}": {
                                "type": "map",
                                "description": "Groovy Map containing sample information. e.g. `[ id:'sample1' ]`",
                            }
                        }
                        inputs_array.append([meta_entry, channel_entry])
                    else:
                        inputs_array.append(channel_entry)
                meta_yml["input"] = ruamel.yaml.comments.CommentedSeq(inputs_array)
                meta_yml["input"].yaml_set_start_comment(
                    "# TODO nf-core: Update the information obtained from bio.tools and make sure that it is correct"
                )
            elif not self.has_meta:
                meta_yml["input"] = [
                    {
                        "bam": {
                            "type": "file",
                            "description": "Sorted BAM/CRAM/SAM file",
                            "pattern": "*.{bam,cram,sam}",
                            "ontologies": [
                                ruamel.yaml.comments.CommentedMap({"edam": "http://edamontology.org/format_2572"}),
                                ruamel.yaml.comments.CommentedMap({"edam": "http://edamontology.org/format_2573"}),
                                ruamel.yaml.comments.CommentedMap({"edam": "http://edamontology.org/format_3462"}),
                            ],
                        }
                    }
                ]
                meta_yml["input"][0]["bam"]["ontologies"][0].yaml_add_eol_comment("BAM", "edam")
                meta_yml["input"][0]["bam"]["ontologies"][1].yaml_add_eol_comment("CRAM", "edam")
                meta_yml["input"][0]["bam"]["ontologies"][2].yaml_add_eol_comment("SAM", "edam")

            if hasattr(self, "outputs") and len(self.outputs) > 0:
                outputs_dict: dict[str, list | dict] = {}
                for i, (output_name, ontologies) in enumerate(self.outputs.items()):
                    channel_contents: list[list[dict] | dict] = []
                    if self.has_meta:
                        channel_contents.append(
                            [
                                {
                                    "meta": {
                                        "type": "map",
                                        "description": "Groovy Map containing sample information. e.g. `[ id:'sample1' ]`",
                                    }
                                }
                            ]
                        )
                    pattern = f"*.{{{','.join(ontologies[2])}}}"
                    file_entry: dict[str, dict] = {
                        pattern: {
                            "type": "file",
                            "description": f"{output_name} file",
                            "pattern": pattern,
                            "ontologies": [
                                ruamel.yaml.comments.CommentedMap({"edam": f"{ont_id}"}) for ont_id in ontologies[0]
                            ],
                        }
                    }
                    for j, ont_desc in enumerate(ontologies[1]):
                        file_entry[pattern]["ontologies"][j].yaml_add_eol_comment(ont_desc, "edam")
                    if self.has_meta:
                        if isinstance(channel_contents[0], list):  # for mypy
                            channel_contents[0].append(file_entry)
                    else:
                        channel_contents.append(file_entry)
                    outputs_dict[output_name] = channel_contents
                outputs_dict.update(versions)
                meta_yml["output"] = ruamel.yaml.comments.CommentedMap(outputs_dict)
                meta_yml["output"].yaml_set_start_comment(
                    "# TODO nf-core: Update the information obtained from bio.tools and make sure that it is correct"
                )
            elif not self.has_meta:
                meta_yml["output"] = {
                    "bam": [
                        {
                            "*.bam": {
                                "type": "file",
                                "description": "Sorted BAM/CRAM/SAM file",
                                "pattern": "*.{bam,cram,sam}",
                                "ontologies": [
                                    ruamel.yaml.comments.CommentedMap({"edam": "http://edamontology.org/format_2572"}),
                                    ruamel.yaml.comments.CommentedMap({"edam": "http://edamontology.org/format_2573"}),
                                    ruamel.yaml.comments.CommentedMap({"edam": "http://edamontology.org/format_3462"}),
                                ],
                            }
                        }
                    ]
                }
                meta_yml["output"]["bam"][0]["*.bam"]["ontologies"][0].yaml_add_eol_comment("BAM", "edam")
                meta_yml["output"]["bam"][0]["*.bam"]["ontologies"][1].yaml_add_eol_comment("CRAM", "edam")
                meta_yml["output"]["bam"][0]["*.bam"]["ontologies"][2].yaml_add_eol_comment("SAM", "edam")
                meta_yml["output"].update(versions)

            meta_yml["topics"] = versions_topic

        else:
            input_entry: list[dict] = [
                {"input": {"type": "file", "description": "", "pattern": "", "ontologies": [{"edam": ""}]}}
            ]
            output_entry: list[dict] = [
                {"*": {"type": "file", "description": "", "pattern": "", "ontologies": [{"edam": ""}]}}
            ]
            if self.has_meta:
                empty_meta_entry: list[dict] = [
                    {
                        "meta": {
                            "type": "map",
                            "description": "Groovy Map containing sample information. e.g. `[ id:'sample1' ]`",
                        }
                    }
                ]
                meta_yml["input"] = [empty_meta_entry + input_entry]
                meta_yml["output"] = {"output": [empty_meta_entry + output_entry]}
            else:
                meta_yml["input"] = input_entry
                meta_yml["output"] = {"output": output_entry}
            meta_yml["output"].update(versions)
            meta_yml["topics"] = versions_topic

        with open(self.file_paths["meta.yml"], "w") as fh:
            yaml.dump(meta_yml, fh)
