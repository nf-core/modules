"""
The NFCoreComponent class holds information and utility functions for a single module or subworkflow
"""

import logging
import re
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)


class NFCoreComponent:
    """
    A class to hold the information about a nf-core module or subworkflow.
    Includes functionality for linting.
    """

    def __init__(
        self,
        component_name: str,
        repo_url: str | None,
        component_dir: Path,
        repo_type: str | None,
        base_dir: Path,
        component_type: str,
        remote_component: bool = True,
    ):
        """
        Initialize the object

        Args:
            component_name (str): The name of the module or subworkflow
            repo_url (str): The URL of the repository
            component_dir (Path): The absolute path to the module or subworkflow
            repo_type (str): Either 'pipeline' or 'modules' depending on
                             whether the directory is a pipeline or clone
                             of nf-core/modules.
            base_dir (Path): The absolute path to the pipeline base dir
            component_type (str): Either 'modules' or 'subworkflows'
            remote_component (bool): Whether the module is to be treated as a
                                     nf-core or local component
        """
        self.component_type = component_type
        self.component_name = component_name
        self.repo_url = repo_url
        self.component_dir = component_dir
        self.repo_type = repo_type
        self.base_dir = base_dir
        self.passed: list[tuple[str, str, str, Path]] = []
        self.warned: list[tuple[str, str, str, Path]] = []
        self.failed: list[tuple[str, str, str, Path]] = []
        self.inputs: list[list[dict[str, dict[str, str]]]] = []
        self.outputs: list[str] = []
        self.topics: dict[str, list[dict[str, dict] | list[dict[str, dict[str, str]]]]]
        self.has_meta: bool = False
        self.git_sha: str | None = None
        self.is_patched: bool = False
        self.branch: str | None = None
        self.workflow_name: str | None = None

        if remote_component:
            # Initialize the important files
            self.main_nf: Path = Path(self.component_dir, "main.nf")
            self.meta_yml: Path | None = Path(self.component_dir, "meta.yml")
            self.environment_yml: Path | None = Path(self.component_dir, "environment.yml")

            component_list = self.component_name.split("/")

            name_index = len(self.component_dir.parts) - 1 - self.component_dir.parts[::-1].index(component_list[0])
            if len(component_list) != 1 and component_list[0] == component_list[1]:
                # Handle cases where the subtool has the same name as the tool
                name_index -= 1

            repo_dir = self.component_dir.parts[:name_index][-1]

            self.org = repo_dir
            self.nftest_testdir: Path | None = Path(self.component_dir, "tests")
            self.nftest_main_nf: Path | None = Path(self.nftest_testdir, "main.nf.test")

            if self.repo_type == "pipeline":
                patch_fn = f"{self.component_name.replace('/', '-')}.diff"
                patch_path = Path(self.component_dir, patch_fn)
                if patch_path.exists():
                    self.is_patched = True
                    self.patch_path = patch_path
        else:
            # The main file is just the local module
            if self.component_dir.is_dir():
                self.main_nf = Path(self.component_dir, "main.nf")
                self.component_name = self.component_dir.stem
                # These attributes are only required by nf-core modules
                # so just set them to None if they don't exist
                self.meta_yml = p if (p := Path(self.component_dir, "meta.yml")).exists() else None
                self.environment_yml = p if (p := Path(self.component_dir, "environment.yml")).exists() else None
                self.nftest_testdir = p if (p := Path(self.component_dir, "tests")).exists() else None
                if self.nftest_testdir is not None:
                    self.nftest_main_nf = p if (p := Path(self.nftest_testdir, "main.nf.test")).exists() else None
            else:
                self.main_nf = self.component_dir
                self.component_dir = self.component_dir.parent
                self.meta_yml = None
                self.environment_yml = None
                self.nftest_testdir = None
                self.nftest_main_nf = None

        self.process_name: str = self._get_process_name()

    def __repr__(self) -> str:
        return f"<NFCoreComponent {self.component_name} {self.component_dir} {self.repo_url}>"

    def _get_main_nf_tags(self, test_main_nf: Path | str):
        """Collect all tags from the main.nf.test file."""
        tags = []
        with open(test_main_nf) as fh:
            for line in fh:
                if line.strip().startswith("tag"):
                    tags.append(line.strip().split()[1].strip('"'))
        return tags

    def _get_included_components(self, main_nf: Path | str):
        """Collect all included components from the main.nf file."""
        included_components = []
        with open(main_nf) as fh:
            for line in fh:
                if line.strip().startswith("include"):
                    # get tool/subtool or subworkflow name from include statement, can be in the form
                    #'../../../modules/nf-core/hisat2/align/main'
                    #'../bam_sort_stats_samtools/main'
                    #'../subworkflows/nf-core/bam_sort_stats_samtools/main'
                    #'plugin/nf-validation'
                    component = line.strip().split()[-1].split(self.org)[-1].split("main")[0].strip("/")
                    component = component.replace("'../", "subworkflows/")
                    component = component.replace("'", "")
                    included_components.append(component)
        return included_components

    def _get_included_components_in_chained_tests(self, main_nf_test: Path | str):
        """Collect all included components from the main.nf file."""
        included_components = []
        with open(main_nf_test) as fh:
            for line in fh:
                if line.strip().startswith("script"):
                    # get tool/subtool or subworkflow name from script statement, can be:
                    # if the component is a module TOOL/SUBTOOL:
                    # '../../SUBTOOL/main.nf'
                    # '../../../TOOL/SUBTOOL/main.nf'
                    # '../../../TOOL/main.nf'
                    # if the component is a module TOOL:
                    # '../../TOOL/main.nf'
                    # '../../TOOL/SUBTOOL/main.nf'
                    # if the component uses full paths or is a subworkflow:
                    # '(../../)modules/nf-core/TOOL/(SUBTOOL/)main.nf'
                    # '(../../)subworkflows/nf-core/TOOL/(SUBTOOL/)main.nf'
                    # the line which uses the current component script:
                    # '../main.nf'
                    component = (
                        line.strip()
                        .split("../")[-1]
                        .split(self.org)[-1]
                        .split("main.nf")[0]
                        .strip("'")
                        .strip('"')
                        .strip("/")
                    )
                    if (
                        "/" in self.component_name
                        and "/" not in component
                        and line.count("../") == 2
                        and self.org not in line
                        and component != ""
                    ):
                        # Add the current component name "TOOL" to the tag
                        component = f"{self.component_name.split('/')[0]}/{component}"
                    if "subworkflows" in line:
                        # Add the subworkflows prefix to the tag
                        component = f"subworkflows/{component}"
                    if component != "":
                        included_components.append(component)
        return included_components

    def _get_process_name(self):
        with open(self.main_nf) as fh:
            for line in fh:
                if re.search(r"^\s*process\s*\w*\s*{", line):
                    return re.search(r"^\s*process\s*(\w*)\s*{.*", line).group(1) or ""
        return ""

    def get_inputs_from_main_nf(self) -> None:
        """Collect all inputs from the main.nf file."""
        inputs: Any = []  # Can be 'list[list[dict[str, dict[str, str]]]]' or 'list[str]'
        with open(self.main_nf) as f:
            data = f.read()
        if self.component_type == "modules":
            # get input values from main.nf after "input:", which can be formatted as tuple val(foo) path(bar) or val foo or val bar or path bar or path foo
            # regex matches:
            # val(foo)
            # path(bar)
            # val foo
            # val bar
            # path bar
            # path foo
            # don't match anything inside comments or after "output:"
            if "input:" not in data:
                log.debug(f"Could not find any inputs in {self.main_nf}")
                return
            input_data = data.split("input:")[1].split("output:")[0]
            for line in input_data.split("\n"):
                channel_elements: Any = []
                line = line.split("//")[0]  # remove any trailing comments
                regex = r"\b(val|path)\b\s*(\(([^)]+)\)|\s*([^)\s,]+))"
                matches = re.finditer(regex, line)
                for _, match in enumerate(matches, start=1):
                    input_val = None
                    if match.group(3):
                        input_val = match.group(3).split(",")[0]  # handle `files, stageAs: "inputs/*"` cases
                    elif match.group(4):
                        input_val = match.group(4).split(",")[0]  # handle `files, stageAs: "inputs/*"` cases
                    if input_val:
                        input_val = re.split(r',(?=(?:[^\'"]*[\'"][^\'"]*[\'"])*[^\'"]*$)', input_val)[
                            0
                        ]  # Takes only first part, avoid commas in quotes
                        input_val = input_val.strip().strip("'").strip('"')  # remove quotes and whitespaces
                        channel_elements.append({input_val: {}})
                if len(channel_elements) == 1:
                    inputs.append(channel_elements[0])
                elif len(channel_elements) > 1:
                    inputs.append(channel_elements)
            log.debug(f"Found {len(inputs)} inputs in {self.main_nf}")
            log.debug(f"Inputs: {inputs}")
            self.inputs = inputs
        elif self.component_type == "subworkflows":
            # get input values from main.nf after "take:"
            if "take:" not in data:
                log.debug(f"Could not find any inputs in {self.main_nf}")
                return
            # get all lines between "take" and "main" or "emit"
            input_data = data.split("take:")[1].split("main:")[0].split("emit:")[0]
            for line in input_data.split("\n"):
                try:
                    inputs.append(line.split()[0])
                except IndexError:
                    pass  # Empty lines
            log.debug(f"Found {len(inputs)} inputs in {self.main_nf}")
            self.inputs = inputs

    def get_outputs_from_main_nf(self):
        with open(self.main_nf) as f:
            data = f.read()
        if self.component_type == "modules":
            outputs = {}
            # get output values from main.nf after "output:". the names are always after "emit:"
            if "output:" not in data:
                log.debug(f"Could not find any outputs in {self.main_nf}")
                return outputs
            output_data = data.split("output:")[1].split("when:")[0]
            log.debug(f"Found output_data: {output_data}")
            regex_emit = r"emit:\s*([^)\s,]+)"
            regex_elements = r"\b(val|path|env|stdout|eval)\b\s*(\(([^)]+)\)|\s*([^)\s,]+))"
            for line in output_data.split("\n"):
                match_emit = re.search(regex_emit, line)
                matches_elements = re.finditer(regex_elements, line)
                if not match_emit:
                    continue
                channel_elements = []
                outputs[match_emit.group(1)] = []
                for _, match_element in enumerate(matches_elements, start=1):
                    output_val = None
                    if match_element.group(3):
                        output_val = match_element.group(3)
                    elif match_element.group(4):
                        output_val = match_element.group(4)
                    if output_val:
                        output_val = re.split(r',(?=(?:[^\'"]*[\'"][^\'"]*[\'"])*[^\'"]*$)', output_val)[
                            0
                        ]  # Takes only first part, avoid commas in quotes
                        output_val = output_val.strip().strip("'").strip('"')  # remove quotes and whitespaces
                        channel_elements.append({output_val: {}})
                if len(channel_elements) == 1:
                    outputs[match_emit.group(1)].append(channel_elements[0])
                elif len(channel_elements) > 1:
                    outputs[match_emit.group(1)].append(channel_elements)
            log.debug(f"Found {len(list(outputs.keys()))} outputs in {self.main_nf}")
            log.debug(f"Outputs: {outputs}")
            self.outputs = outputs
        elif self.component_type == "subworkflows":
            outputs = []
            # get output values from main.nf after "emit:". Can be named outputs or not.
            if "emit:" not in data:
                log.debug(f"Could not find any outputs in {self.main_nf}")
                return outputs
            output_data = data.split("emit:")[1].split("}")[0]
            for line in output_data.split("\n"):
                try:
                    outputs.append(line.split("=")[0].split()[0])
                except IndexError:
                    # Empty lines
                    pass
            log.debug(f"Found {len(outputs)} outputs in {self.main_nf}")
            self.outputs = outputs

    def get_topics_from_main_nf(self) -> None:
        with open(self.main_nf) as f:
            data = f.read()
        if self.component_type == "modules":
            topics: dict[str, list[dict[str, dict] | list[dict[str, dict[str, str]]]]] = {}
            # get topic name from main.nf after "output:". the names are always after "topic:"
            if "output:" not in data:
                log.debug(f"Could not find any outputs in {self.main_nf}")
                self.topics = topics
                return
            output_data = data.split("output:")[1].split("when:")[0]
            log.debug(f"Output data: {output_data}")
            regex_topic = r"topic:\s*([^)\s,]+)"
            regex_elements = r"\b(val|path|env|stdout|eval)\b\s*(\(([^)]+)\)|\s*([^)\s,]+))"
            for line in output_data.split("\n"):
                match_topic = re.search(regex_topic, line)
                matches_elements = re.finditer(regex_elements, line)
                if not match_topic:
                    continue
                channel_elements: list[dict[str, dict]] = []
                topic_name = match_topic.group(1)
                if topic_name in topics:
                    continue
                topics[match_topic.group(1)] = []
                for _, match_element in enumerate(matches_elements, start=1):
                    topic_val = None
                    if match_element.group(3):
                        topic_val = match_element.group(3)
                    elif match_element.group(4):
                        topic_val = match_element.group(4)
                    if topic_val:
                        topic_val = re.split(r',(?=(?:[^\'"]*[\'"][^\'"]*[\'"])*[^\'"]*$)', topic_val)[
                            0
                        ]  # Takes only first part, avoid commas in quotes
                        topic_val = topic_val.strip().strip("'").strip('"')  # remove quotes and whitespaces
                        channel_elements.append({topic_val: {}})
                if len(channel_elements) == 1:
                    topics[match_topic.group(1)].append(channel_elements[0])
                elif len(channel_elements) > 1:
                    topics[match_topic.group(1)].append(channel_elements)
            log.debug(f"Found {len(list(topics.keys()))} topics in {self.main_nf}")
            log.debug(f"Topics: {topics}")
            self.topics = topics
