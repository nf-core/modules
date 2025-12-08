"""
Code for linting subworkflows in the nf-core/subworkflows repository and
in nf-core pipelines

Command:
nf-core subworkflows lint
"""

import logging
import os

import questionary
import rich
import ruamel.yaml

import nf_core.utils
from nf_core.components.lint import ComponentLint, LintExceptionError, LintResult
from nf_core.pipelines.lint_utils import console, run_prettier_on_file

log = logging.getLogger(__name__)

# Import lint functions
from .main_nf import main_nf  # type: ignore[misc]
from .meta_yml import meta_yml  # type: ignore[misc]
from .subworkflow_changes import subworkflow_changes  # type: ignore[misc]
from .subworkflow_if_empty_null import subworkflow_if_empty_null  # type: ignore[misc]
from .subworkflow_tests import subworkflow_tests  # type: ignore[misc]
from .subworkflow_todos import subworkflow_todos  # type: ignore[misc]
from .subworkflow_version import subworkflow_version  # type: ignore[misc]


class SubworkflowLint(ComponentLint):
    """
    An object for linting subworkflows either in a clone of the 'nf-core/modules'
    repository or in any nf-core pipeline directory
    """

    main_nf = main_nf
    meta_yml = meta_yml
    subworkflow_changes = subworkflow_changes
    subworkflow_tests = subworkflow_tests
    subworkflow_todos = subworkflow_todos
    subworkflow_if_empty_null = subworkflow_if_empty_null
    subworkflow_version = subworkflow_version

    def __init__(
        self,
        directory,
        fail_warned=False,
        fix=False,
        remote_url=None,
        branch=None,
        no_pull=False,
        registry=None,
        hide_progress=False,
    ):
        super().__init__(
            component_type="subworkflows",
            directory=directory,
            fail_warned=fail_warned,
            fix=fix,
            remote_url=remote_url,
            branch=branch,
            no_pull=no_pull,
            registry=registry,
            hide_progress=hide_progress,
        )

    def lint(
        self,
        subworkflow=None,
        registry="quay.io",
        key=(),
        all_subworkflows=False,
        print_results=True,
        show_passed=False,
        sort_by="test",
        local=False,
    ):
        """
        Lint all or one specific subworkflow

        First gets a list of all local subworkflows (in subworkflows/local/process) and all subworkflows
        installed from nf-core (in subworkflows/nf-core)

        For all nf-core subworkflows, the correct file structure is assured and important
        file content is verified. If directory subject to linting is a clone of 'nf-core/modules',
        the files necessary for testing the subworkflows are also inspected.

        For all local subworkflows, the '.nf' file is checked for some important flags, and warnings
        are issued if untypical content is found.

        :param subworkflow:     A specific subworkflow to lint
        :param print_results:   Whether to print the linting results
        :param show_passed:     Whether passed tests should be shown as well
        :param hide_progress:   Don't show progress bars

        :returns:               A SubworkflowLint object containing information of
                                the passed, warned and failed tests
        """
        # TODO: consider unifying modules and subworkflows lint() function and add it to the ComponentLint class
        # Prompt for subworkflow or all
        if subworkflow is None and not (local or all_subworkflows):
            questions = [
                {
                    "type": "list",
                    "name": "all_subworkflows",
                    "message": "Lint all subworkflows or a single named subworkflow?",
                    "choices": ["All subworkflows", "Named subworkflow"],
                },
                {
                    "type": "autocomplete",
                    "name": "subworkflow_name",
                    "message": "Subworkflow name:",
                    "when": lambda x: x["all_subworkflows"] == "Named subworkflow",
                    "choices": [m.component_name for m in self.all_remote_components],
                },
            ]
            answers = questionary.unsafe_prompt(questions, style=nf_core.utils.nfcore_question_style)
            all_subworkflows = answers["all_subworkflows"] == "All subworkflows"
            subworkflow = answers.get("subworkflow_name")

        # Only lint the given module
        if subworkflow:
            if all_subworkflows:
                raise LintExceptionError("You cannot specify a tool and request all tools to be linted.")
            local_subworkflows = []
            remote_subworkflows = [s for s in self.all_remote_components if s.component_name == subworkflow]
            if len(remote_subworkflows) == 0:
                raise LintExceptionError(f"Could not find the specified subworkflow: '{subworkflow}'")
        else:
            local_subworkflows = self.all_local_components
            remote_subworkflows = self.all_remote_components

        if self.repo_type == "modules":
            log.info(f"Linting modules repo: [magenta]'{self.directory}'")
        else:
            log.info(f"Linting pipeline: [magenta]'{self.directory}'")
        if subworkflow:
            log.info(f"Linting subworkflow: [magenta]'{subworkflow}'")

        # Filter the tests by the key if one is supplied
        if key:
            self.filter_tests_by_key(key)
            log.info("Only running tests: '{}'".format("', '".join(key)))

        # If it is a pipeline, load the lint config file and the modules.json file
        if self.repo_type == "pipeline":
            self.set_up_pipeline_files()

        # Lint local subworkflows
        if local and len(local_subworkflows) > 0:
            self.lint_subworkflows(local_subworkflows, registry=registry, local=True)

        # Lint nf-core subworkflows
        if not local and len(remote_subworkflows) > 0:
            self.lint_subworkflows(remote_subworkflows, registry=registry, local=False)

        if print_results:
            self._print_results(show_passed=show_passed, sort_by=sort_by)
            self.print_summary()

    def lint_subworkflows(self, subworkflows, registry="quay.io", local=False):
        """
        Lint a list of subworkflows

        Args:
            subworkflows ([NFCoreComponent]): A list of subworkflow objects
            registry (str): The container registry to use. Should be quay.io in most situations.
            local (boolean): Whether the list consist of local or nf-core subworkflows
        """
        # TODO: consider unifying modules and subworkflows lint_subworkflows() function and add it to the ComponentLint class
        progress_bar = rich.progress.Progress(
            "[bold blue]{task.description}",
            rich.progress.BarColumn(bar_width=None),
            "[magenta]{task.completed} of {task.total}[reset] » [bold yellow]{task.fields[test_name]}",
            transient=True,
            console=console,
            disable=self.hide_progress or os.environ.get("HIDE_PROGRESS", None) is not None,
        )
        with progress_bar:
            lint_progress = progress_bar.add_task(
                f"Linting {'local' if local else 'nf-core'} subworkflows",
                total=len(subworkflows),
                test_name=subworkflows[0].component_name,
            )

            for swf in subworkflows:
                progress_bar.update(lint_progress, advance=1, test_name=swf.component_name)
                self.lint_subworkflow(swf, progress_bar, registry=registry, local=local)

    def lint_subworkflow(self, swf, progress_bar, registry, local=False):
        """
        Perform linting on one subworkflow

        If the subworkflow is a local subworkflow we only check the `main.nf` file,
        and issue warnings instead of failures.

        If the subworkflow is a nf-core subworkflow we check for existence of the files
        - main.nf
        - meta.yml
        And verify that their content conform to the nf-core standards.

        If the linting is run for subworkflows in the central nf-core/modules repo
        (repo_type==modules), files that are relevant for subworkflow testing are
        also examined
        """
        # TODO: consider unifying modules and subworkflows lint_subworkflow() function and add it to the ComponentLint class
        # Only check the main script in case of a local subworkflow
        if local:
            self.main_nf(swf)
            self.meta_yml(swf, allow_missing=True)
            self.subworkflow_todos(swf)
            self.passed += [LintResult(swf, *s) for s in swf.passed]
            warned = [LintResult(swf, *m) for m in (swf.warned + swf.failed)]
            if not self.fail_warned:
                self.warned += warned
            else:
                self.failed += warned

        # Otherwise run all the lint tests
        else:
            # Update meta.yml file if requested
            if self.fix:
                self.update_meta_yml_file(swf)

            if self.repo_type == "pipeline" and self.modules_json:
                # Set correct sha
                version = self.modules_json.get_subworkflow_version(swf.component_name, swf.repo_url, swf.org)
                swf.git_sha = version

            for test_name in self.lint_tests:
                getattr(self, test_name)(swf)

            self.passed += [LintResult(swf, *s) for s in swf.passed]
            warned = [LintResult(swf, *s) for s in swf.warned]
            if not self.fail_warned:
                self.warned += warned
            else:
                self.failed += warned

            self.failed += [LintResult(swf, *s) for s in swf.failed]

    def update_meta_yml_file(self, swf):
        """
        Update the meta.yml file with the correct inputs and outputs
        """
        yaml = ruamel.yaml.YAML()
        yaml.preserve_quotes = True
        yaml.indent(mapping=2, sequence=2, offset=0)

        # Read meta.yml
        with open(swf.meta_yml) as fh:
            meta_yaml = yaml.load(fh)
        meta_yaml_corrected = meta_yaml.copy()
        # Obtain inputs and outputs from main.nf
        swf.get_inputs_from_main_nf()
        swf.get_outputs_from_main_nf()

        # Compare inputs and add them if missing
        if "input" in meta_yaml:
            # Delete inputs from meta.yml which are not present in main.nf
            meta_yaml_corrected["input"] = [
                input for input in meta_yaml["input"] if list(input.keys())[0] in swf.inputs
            ]
            # Obtain inputs from main.nf missing in meta.yml
            inputs_correct = [
                list(input.keys())[0] for input in meta_yaml_corrected["input"] if list(input.keys())[0] in swf.inputs
            ]
            inputs_missing = [input for input in swf.inputs if input not in inputs_correct]
            # Add missing inputs to meta.yml
            for missing_input in inputs_missing:
                meta_yaml_corrected["input"].append({missing_input: {"description": ""}})

        if "output" in meta_yaml:
            # Delete outputs from meta.yml which are not present in main.nf
            meta_yaml_corrected["output"] = [
                output for output in meta_yaml["output"] if list(output.keys())[0] in swf.outputs
            ]
            # Obtain output from main.nf missing in meta.yml
            outputs_correct = [
                list(output.keys())[0]
                for output in meta_yaml_corrected["output"]
                if list(output.keys())[0] in swf.outputs
            ]
            outputs_missing = [output for output in swf.outputs if output not in outputs_correct]
            # Add missing inputs to meta.yml
            for missing_output in outputs_missing:
                meta_yaml_corrected["output"].append({missing_output: {"description": ""}})

        # Write corrected meta.yml to file
        with open(swf.meta_yml, "w") as fh:
            log.info(f"Updating {swf.meta_yml}")
            yaml.dump(meta_yaml_corrected, fh)
            run_prettier_on_file(fh.name)
