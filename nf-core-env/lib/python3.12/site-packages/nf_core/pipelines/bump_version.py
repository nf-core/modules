"""Bumps the version number in all appropriate files for
a nf-core pipeline.
"""

import logging
import re
from pathlib import Path

import rich.console
from ruamel.yaml import YAML

import nf_core.utils
from nf_core.pipelines.rocrate import ROCrate
from nf_core.utils import Pipeline

log = logging.getLogger(__name__)
stderr = rich.console.Console(stderr=True, force_terminal=nf_core.utils.rich_force_colors())


def bump_pipeline_version(pipeline_obj: Pipeline, new_version: str) -> None:
    """Bumps a pipeline version number.

    Args:
        pipeline_obj (nf_core.utils.Pipeline): A `Pipeline` object that holds information
            about the pipeline contents and build files.
        new_version (str): The new version tag for the pipeline. Semantic versioning only.
    """

    # Collect the old and new version numbers
    current_version = pipeline_obj.nf_config.get("manifest.version", "").strip(" '\"")
    if new_version.startswith("v"):
        log.warning("Stripping leading 'v' from new version number")
        new_version = new_version[1:]
    if not current_version:
        raise UserWarning("Could not find config variable 'manifest.version'")
    if current_version == new_version:
        raise UserWarning(f"Current version is already: {current_version}")
    log.info(f"Changing version number from '{current_version}' to '{new_version}'")

    # nextflow.config - workflow manifest version
    update_file_version(
        "nextflow.config",
        pipeline_obj,
        [
            (
                rf"(version\s*=\s*['\"]){re.escape(current_version)}(['\"])",
                rf"\g<1>{new_version}\g<2>",
            )
        ],
    )
    # multiqc_config.yaml
    multiqc_new_version = "dev" if "dev" in new_version else new_version
    multiqc_current_version = "dev" if "dev" in current_version else current_version
    if multiqc_current_version != "dev" and multiqc_new_version != "dev":
        update_file_version(
            Path("assets", "multiqc_config.yml"),
            pipeline_obj,
            [
                (
                    f"/releases/tag/{current_version}",
                    f"/releases/tag/{new_version}",
                )
            ],
            yaml_key=["report_comment"],
        )
    if multiqc_current_version != "dev" and multiqc_new_version == "dev":
        update_file_version(
            Path("assets", "multiqc_config.yml"),
            pipeline_obj,
            [
                (
                    f"/releases/tag/{current_version}",
                    "/tree/dev",
                )
            ],
            yaml_key=["report_comment"],
        )
    if multiqc_current_version == "dev" and multiqc_new_version != "dev":
        update_file_version(
            Path("assets", "multiqc_config.yml"),
            pipeline_obj,
            [
                (
                    "/tree/dev",
                    f"/releases/tag/{multiqc_new_version}",
                )
            ],
            yaml_key=["report_comment"],
        )
    update_file_version(
        Path("assets", "multiqc_config.yml"),
        pipeline_obj,
        [
            (
                f"/{multiqc_current_version}/",
                f"/{multiqc_new_version}/",
            ),
        ],
        yaml_key=["report_comment"],
    )
    # nf-test snap files
    pipeline_name = pipeline_obj.nf_config.get("manifest.name", "").strip(" '\"")
    snap_files = [f.relative_to(pipeline_obj.wf_path) for f in Path(pipeline_obj.wf_path).glob("tests/pipeline/*.snap")]
    for snap_file in snap_files:
        update_file_version(
            snap_file,
            pipeline_obj,
            [
                (
                    f"{pipeline_name}={current_version}",
                    f"{pipeline_name}={new_version}",
                )
            ],
            required=False,
        )
    # .nf-core.yml - pipeline version
    # update entry: version: 1.0.0dev, but not `nf_core_version`, or `bump_version`
    update_file_version(
        ".nf-core.yml",
        pipeline_obj,
        [
            (
                current_version,
                new_version,
            )
        ],
        required=False,
        yaml_key=["template", "version"],
    )

    # update rocrate if ro-crate is present
    if Path(pipeline_obj.wf_path, "ro-crate-metadata.json").exists():
        ROCrate(pipeline_obj.wf_path).update_rocrate()


def bump_nextflow_version(pipeline_obj: Pipeline, new_version: str) -> None:
    """Bumps the required Nextflow version number of a pipeline.

    Args:
        pipeline_obj (nf_core.utils.Pipeline): A `Pipeline` object that holds information
            about the pipeline contents and build files.
        new_version (str): The new version tag for the required Nextflow version.
    """

    # Collect the old and new version numbers - strip leading non-numeric characters (>=)
    current_version = pipeline_obj.nf_config.get("manifest.nextflowVersion", "").strip(" '\"")
    current_version = re.sub(r"^[^0-9\.]*", "", current_version)
    new_version = re.sub(r"^[^0-9\.]*", "", new_version)
    if not current_version:
        raise UserWarning("Could not find config variable 'manifest.nextflowVersion'")
    log.info(f"Changing Nextlow version number from '{current_version}' to '{new_version}'")

    # nextflow.config - manifest minimum nextflowVersion
    update_file_version(
        "nextflow.config",
        pipeline_obj,
        [
            (
                rf"(nextflowVersion\s*=\s*[\'\"]?!>=\s*)({re.escape(current_version)})([\'\"]?)",
                rf"\g<1>{new_version}\g<3>",
            )
        ],
    )

    # .github/workflows/nf-test.yml - Nextflow version matrix
    update_file_version(
        Path(".github", "workflows", "nf-test.yml"),
        pipeline_obj,
        [
            (
                # example:
                # NXF_VER:
                #   - "20.04.0"
                current_version,
                new_version,
            )
        ],
        yaml_key=["jobs", "nf-test", "strategy", "matrix", "NXF_VER"],
    )

    # README.md - Nextflow version badge
    update_file_version(
        "README.md",
        pipeline_obj,
        [
            (
                rf"nextflow%20DSL2-%E2%89%A5{re.escape(current_version)}-23aa62.svg",
                f"nextflow%20DSL2-%E2%89%A5{new_version}-23aa62.svg",
            ),
            (f"version-%E2%89%A5{re.escape(current_version)}-green", f"version-%E2%89%A5{new_version}-green"),
        ],
        False,
    )


def update_file_version(
    filename: str | Path,
    pipeline_obj: Pipeline,
    patterns: list[tuple[str, str]],
    required: bool = True,
    yaml_key: list[str] | None = None,
) -> None:
    """
    Updates a file with a new version number.

    Args:
        filename (str): The name of the file to update.
        pipeline_obj (nf_core.utils.Pipeline): A `Pipeline` object that holds information
            about the pipeline contents.
        patterns (List[Tuple[str, str]]): A list of tuples containing the regex patterns to
            match and the replacement strings.
        required (bool, optional): Whether the file is required to exist. Defaults to `True`.
        yaml_key (List[str] | None, optional): The YAML key to update. Defaults to `None`.
    """
    fn: Path = pipeline_obj._fp(filename)

    if not fn.exists():
        log.warning(f"File not found: '{fn}'")
        return

    if yaml_key:
        update_yaml_file(fn, patterns, yaml_key, required)
    else:
        update_text_file(fn, patterns, required)


def update_yaml_file(fn: Path, patterns: list[tuple[str, str]], yaml_key: list[str], required: bool):
    """
    Updates a YAML file with a new version number.

    Args:
        fn (Path): The name of the file to update.
        patterns (List[Tuple[str, str]]): A list of tuples containing the regex patterns to
            match and the replacement strings.
        yaml_key (List[str]): The YAML key to update.
        required (bool): Whether the file is required to exist.
    """
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.width = 4096  # Prevent line wrapping
    with open(fn) as file:
        yaml_content = yaml.load(file)

    try:
        target = yaml_content
        for key in yaml_key[:-1]:
            target = target[key]

        last_key = yaml_key[-1]
        current_value = target[last_key]

        new_value = current_value
        for pattern, replacement in patterns:
            # check if current value is list
            if isinstance(current_value, list):
                new_value = [re.sub(pattern, replacement, item) for item in current_value]
            else:
                new_value = re.sub(pattern, replacement, current_value)

        if new_value != current_value:
            target[last_key] = new_value
            with open(fn, "w") as file:
                yaml.indent(mapping=2, sequence=4, offset=2)  # from https://stackoverflow.com/a/44389139/1696643
                yaml.dump(yaml_content, file)
            log.info(f"Updated version in YAML file '{fn}'")
            log_change(str(current_value), str(new_value))
    except KeyError as e:
        handle_error(f"Could not find key {e} in the YAML structure of {fn}", required)


def update_text_file(fn: Path, patterns: list[tuple[str, str]], required: bool):
    """
    Updates a text file with a new version number.

    Args:
        fn (Path): The name of the file to update.
        patterns (List[Tuple[str, str]]): A list of tuples containing the regex patterns to
            match and the replacement strings.
        required (bool): Whether the file is required to exist.
    """
    with open(fn) as file:
        content = file.read()

    updated = False
    for pattern, replacement in patterns:
        new_content, count = re.subn(pattern, replacement, content)
        if count > 0:
            log_change(content, new_content)
            content = new_content
            updated = True
            log.info(f"Updated version in '{fn}'")
            log.debug(f"Replaced pattern '{pattern}' with '{replacement}' {count} times")
        else:
            handle_error(f"Could not find version number in {fn}: `{pattern}`", required)

    if updated:
        with open(fn, "w") as file:
            file.write(content)


def handle_error(message: str, required: bool):
    if required:
        raise ValueError(message)
    else:
        log.info(message)


def log_change(old_content: str, new_content: str):
    old_lines = old_content.splitlines()
    new_lines = new_content.splitlines()

    for old_line, new_line in zip(old_lines, new_lines):
        if old_line != new_line:
            stderr.print(f"          [red] - {old_line.strip()}", highlight=False)
            stderr.print(f"          [green] + {new_line.strip()}", highlight=False)

    stderr.print("\n")
