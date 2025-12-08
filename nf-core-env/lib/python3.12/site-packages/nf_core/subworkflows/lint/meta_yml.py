import json
import logging
from pathlib import Path

import jsonschema.validators
import ruamel.yaml

import nf_core.components.components_utils
from nf_core.components.lint import LintExceptionError

log = logging.getLogger(__name__)


def meta_yml(subworkflow_lint_object, subworkflow, allow_missing: bool = False):
    """
    Lint a ``meta.yml`` file

    The lint test checks that the subworkflow has
    a ``meta.yml`` file and that it follows the
    JSON schema defined in the ``subworkflows/yaml-schema.json``
    file in the nf-core/modules repository.

    In addition it checks that the subworkflow name
    and subworkflow input is consistent between the
    ``meta.yml`` and the ``main.nf``.

    Checks that all input and output channels are specified in ``meta.yml``.
    Checks that all included components in ``main.nf`` are specified in ``meta.yml``.

    """
    # Read the meta.yml file
    if subworkflow.meta_yml is None:
        if allow_missing:
            subworkflow.warned.append(
                (
                    "meta_yml",
                    "meta_yml_exists",
                    "Subworkflow `meta.yml` does not exist",
                    Path(subworkflow.component_dir, "meta.yml"),
                )
            )
            return
        raise LintExceptionError("Subworkflow does not have a `meta.yml` file")

    try:
        with open(subworkflow.meta_yml) as fh:
            yaml = ruamel.yaml.YAML(typ="safe")
            meta_yaml = yaml.load(fh)
        subworkflow.passed.append(
            ("meta_yml", "meta_yml_exists", "Subworkflow `meta.yml` exists", subworkflow.meta_yml)
        )
    except FileNotFoundError:
        subworkflow.failed.append(
            ("meta_yml", "meta_yml_exists", "Subworkflow `meta.yml` does not exist", subworkflow.meta_yml)
        )
        return

    # Confirm that the meta.yml file is valid according to the JSON schema
    valid_meta_yml = True
    try:
        with open(Path(subworkflow_lint_object.modules_repo.local_repo_dir, "subworkflows/yaml-schema.json")) as fh:
            schema = json.load(fh)
        jsonschema.validators.validate(instance=meta_yaml, schema=schema)
        subworkflow.passed.append(
            ("meta_yml", "meta_yml_valid", "Subworkflow `meta.yml` is valid", subworkflow.meta_yml)
        )
    except jsonschema.exceptions.ValidationError as e:
        valid_meta_yml = False
        hint = ""
        if len(e.path) > 0:
            hint = f"\nCheck the entry for `{e.path[0]}`."
        if e.message.startswith("None is not of type 'object'") and len(e.path) > 2:
            hint = f"\nCheck that the child entries of {str(e.path[0]) + '.' + str(e.path[2])} are indented correctly."
        subworkflow.failed.append(
            (
                "meta_yml",
                "meta_yml_valid",
                f"The `meta.yml` of the subworkflow {subworkflow.component_name} is not valid: {e.message}.{hint}",
                subworkflow.meta_yml,
            )
        )
        return

    # Confirm that all input and output channels are specified
    if valid_meta_yml:
        if "input" in meta_yaml:
            meta_input = [list(x.keys())[0] for x in meta_yaml["input"]]
            for input in subworkflow.inputs:
                if input in meta_input:
                    subworkflow.passed.append(("meta_yml", "meta_input", f"`{input}` specified", subworkflow.meta_yml))
                else:
                    subworkflow.failed.append(
                        ("meta_yml", "meta_input", f"`{input}` missing in `meta.yml`", subworkflow.meta_yml)
                    )
        else:
            log.debug(f"No inputs specified in subworkflow `main.nf`: {subworkflow.component_name}")

        if "output" in meta_yaml:
            meta_output = [list(x.keys())[0] for x in meta_yaml["output"]]
            for output in subworkflow.outputs:
                if output in meta_output:
                    subworkflow.passed.append(
                        ("meta_yml", "meta_output", f"`{output}` specified", subworkflow.meta_yml)
                    )
                else:
                    subworkflow.failed.append(
                        ("meta_yml", "meta_output", f"`{output}` missing in `meta.yml`", subworkflow.meta_yml)
                    )
        else:
            log.debug(f"No outputs specified in subworkflow `main.nf`: {subworkflow.component_name}")

        # confirm that the name matches the process name in main.nf
        if meta_yaml["name"].upper() == subworkflow.workflow_name:
            subworkflow.passed.append(
                ("meta_yml", "meta_name", "Correct name specified in `meta.yml`", subworkflow.meta_yml)
            )
        else:
            subworkflow.failed.append(
                (
                    "meta_yml",
                    "meta_name",
                    f"Conflicting workflow name between meta.yml (`{meta_yaml['name']}`) and main.nf (`{subworkflow.workflow_name}`)",
                    subworkflow.meta_yml,
                )
            )

        # confirm that all included components in ``main.nf`` are specified in ``meta.yml``
        included_components_ = nf_core.components.components_utils.get_components_to_install(subworkflow.component_dir)
        included_components = included_components_[0] + included_components_[1]
        # join included modules and included subworkflows in a single list
        included_components_names = [component["name"] for component in included_components]
        if "components" in meta_yaml:
            meta_components = [x if isinstance(x, str) else list(x)[0] for x in meta_yaml["components"]]
            for component in set(included_components_names):
                if component in meta_components:
                    subworkflow.passed.append(
                        (
                            "meta_yml",
                            "meta_include",
                            f"Included module/subworkflow `{component}` specified in `meta.yml`",
                            subworkflow.meta_yml,
                        )
                    )
                else:
                    subworkflow.failed.append(
                        (
                            "meta_yml",
                            "meta_include",
                            f"Included module/subworkflow `{component}` missing in `meta.yml`",
                            subworkflow.meta_yml,
                        )
                    )
        if "modules" in meta_yaml:
            subworkflow.failed.append(
                (
                    "meta_yml",
                    "meta_modules_deprecated",
                    "Deprecated section 'modules' found in `meta.yml`, use 'components' instead",
                    subworkflow.meta_yml,
                )
            )
        else:
            subworkflow.passed.append(
                (
                    "meta_yml",
                    "meta_modules_deprecated",
                    "Deprecated section 'modules' not found in `meta.yml`",
                    subworkflow.meta_yml,
                )
            )
