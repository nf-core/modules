import json
import logging
from pathlib import Path

import ruamel.yaml
from jsonschema import exceptions, validators

from nf_core.components.lint import ComponentLint, LintExceptionError
from nf_core.components.nfcore_component import NFCoreComponent

log = logging.getLogger(__name__)

# Configure ruamel.yaml for proper formatting
yaml = ruamel.yaml.YAML()
yaml.indent(mapping=2, sequence=2, offset=2)


def environment_yml(module_lint_object: ComponentLint, module: NFCoreComponent, allow_missing: bool = False) -> None:
    """
    Lint an ``environment.yml`` file.

    The lint test checks that the ``dependencies`` section
    in the environment.yml file is valid YAML and that it
    is sorted alphabetically.
    """
    env_yml = None
    has_schema_header = False
    lines = []

    # Define the schema lines to be added if missing
    schema_lines = [
        "---\n",
        "# yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json\n",
    ]

    #  load the environment.yml file
    if module.environment_yml is None:
        if allow_missing:
            module.warned.append(
                (
                    "environment_yml",
                    "environment_yml_exists",
                    "Module's `environment.yml` does not exist",
                    Path(module.component_dir, "environment.yml"),
                ),
            )
            return
        raise LintExceptionError("Module does not have an `environment.yml` file")
    try:
        # Read the entire file content to handle headers properly
        with open(module.environment_yml) as fh:
            lines = fh.readlines()

        # Check if the first two lines contain schema configuration
        content_start = 0

        if len(lines) >= 2 and lines[0] == "---\n" and lines[1].startswith("# yaml-language-server: $schema="):
            has_schema_header = True
            content_start = 2

        content = "".join(lines[content_start:])  # Skip schema lines when reading content

        # Parse the YAML content
        env_yml = yaml.load(content)
        if env_yml is None:
            raise ruamel.yaml.scanner.ScannerError("Empty YAML file")

        module.passed.append(
            (
                "environment_yml",
                "environment_yml_exists",
                "Module's `environment.yml` exists",
                module.environment_yml,
            )
        )

    except FileNotFoundError:
        # check if the module's main.nf requires a conda environment
        with open(Path(module.component_dir, "main.nf")) as fh:
            main_nf = fh.read()
            if 'conda "${moduleDir}/environment.yml"' in main_nf:
                module.failed.append(
                    (
                        "environment_yml",
                        "environment_yml_exists",
                        "Module's `environment.yml` does not exist",
                        module.environment_yml,
                    )
                )
            else:
                module.passed.append(
                    (
                        "environment_yml",
                        "environment_yml_exists",
                        "Module's `environment.yml` does not exist, but it is also not included in the main.nf",
                        module.environment_yml,
                    )
                )

    # Confirm that the environment.yml file is valid according to the JSON schema
    if env_yml:
        valid_env_yml = False
        try:
            with open(Path(module_lint_object.modules_repo.local_repo_dir, "modules/environment-schema.json")) as fh:
                schema = json.load(fh)
            validators.validate(instance=env_yml, schema=schema)
            module.passed.append(
                (
                    "environment_yml",
                    "environment_yml_valid",
                    "Module's `environment.yml` is valid",
                    module.environment_yml,
                )
            )
            valid_env_yml = True
        except exceptions.ValidationError as e:
            hint = ""
            if len(e.path) > 0:
                hint = f"\nCheck the entry for `{e.path[0]}`."
            if e.schema and isinstance(e.schema, dict) and "message" in e.schema:
                e.message = e.schema["message"]
            module.failed.append(
                (
                    "environment_yml",
                    "environment_yml_valid",
                    f"The `environment.yml` of the module {module.component_name} is not valid: {e.message}.{hint}",
                    module.environment_yml,
                )
            )

        if valid_env_yml:
            # Sort dependencies if they exist
            if "dependencies" in env_yml:
                dicts = []
                others = []

                for term in env_yml["dependencies"]:
                    if isinstance(term, dict):
                        dicts.append(term)
                    else:
                        others.append(term)

                # Sort non-dict dependencies with special handling for pip
                def sort_key(x):
                    # Convert to string for comparison
                    str_x = str(x)
                    # If it's a pip package (but not pip itself), put it after other conda packages
                    if str_x.startswith("pip=") or str_x == "pip":
                        return (1, str_x)  # pip comes after other conda packages
                    else:
                        return (0, str_x)  # regular conda packages come first

                others.sort(key=sort_key)

                # Sort any lists within dict dependencies
                for dict_term in dicts:
                    for value in dict_term.values():
                        if isinstance(value, list):
                            value.sort(key=str)

                # Sort dict dependencies alphabetically
                dicts.sort(key=str)

                # Combine sorted dependencies
                sorted_deps = others + dicts

                # Check if dependencies are already sorted
                is_sorted = env_yml["dependencies"] == sorted_deps and all(
                    not isinstance(term, dict)
                    or all(not isinstance(value, list) or value == sorted(value, key=str) for value in term.values())
                    for term in env_yml["dependencies"]
                )
            else:
                sorted_deps = None
                is_sorted = True

            if is_sorted:
                module.passed.append(
                    (
                        "environment_yml",
                        "environment_yml_sorted",
                        "The dependencies in the module's `environment.yml` are sorted correctly",
                        module.environment_yml,
                    )
                )
            else:
                log.info(
                    f"Dependencies in {module.component_name}'s environment.yml were not sorted. Sorting them now."
                )

                # Update dependencies if they need sorting
                if sorted_deps is not None:
                    env_yml["dependencies"] = sorted_deps

                # Write back to file with headers
                with open(Path(module.component_dir, "environment.yml"), "w") as fh:
                    # If file had a schema header, check if it's pointing to a different URL
                    if has_schema_header and len(lines) >= 2:
                        existing_schema_line = lines[1]
                        # If the existing schema URL is different, update it
                        if not existing_schema_line.endswith("/modules/master/modules/environment-schema.json\n"):
                            fh.writelines(schema_lines)
                        else:
                            # Keep the existing schema lines
                            fh.writelines(lines[:2])
                    else:
                        # No schema header present, add the default one
                        fh.writelines(schema_lines)
                    # Then dump the sorted YAML
                    yaml.dump(env_yml, fh)

                module.passed.append(
                    (
                        "environment_yml",
                        "environment_yml_sorted",
                        "The dependencies in the module's `environment.yml` have been sorted",
                        module.environment_yml,
                    )
                )
