import re
from pathlib import Path
from textwrap import dedent

import pytest
import yaml


def _get_workflow_names():
    """Get all names of all workflows which have a test.yml in the tests directory.

    To do so, recursively finds all test.yml files and parses their content.
    """
    here = Path(__file__).parent.resolve()
    pytest_workflow_files = here.glob("modules/**/test.yml")
    for f in pytest_workflow_files:
        # test_config = yaml.safe_load(f.read_text())
        test_config = yaml.load(f.read_text(), Loader=yaml.BaseLoader)
        for workflow in test_config:
            # https://github.com/nf-core/modules/pull/1242 - added to cover tests
            # that expect an error and therefore will not generate a versions.yml
            if "exit_code" not in workflow:
                yield workflow["name"]


@pytest.mark.workflow(*_get_workflow_names())
def test_ensure_valid_version_yml(workflow_dir):
    workflow_dir = Path(workflow_dir)
    software_name = workflow_dir.name.split("_")[0].lower()
    try:
        versions_yml_file = workflow_dir / f"output/{software_name}/versions.yml"
        versions_yml = versions_yml_file.read_text()
    except FileNotFoundError:
        raise AssertionError(
            dedent(
                f"""\
                `versions.yml` not found in the output directory.
                Expected path: `{versions_yml_file}`

                This can have multiple reasons:
                * The test-workflow failed before a `versions.yml` could be generated.
                * The workflow name in `test.yml` does not start with the tool name.
                """
            )
        )

    assert (
        "END_VERSIONS" not in versions_yml
    ), "END_VERSIONS detected in versions.yml. This is a sign of an ill-formatted HEREDOC"

    # Raises an exception if yaml is not valid
    versions = yaml.safe_load(versions_yml)
    assert (
        len(versions) == 1
    ), "The top-level of versions.yml must contain exactly one entry: the process name as dict key"
    software_versions = next(iter(versions.values()))
    assert len(software_versions), "There must be at least one version emitted."
    for tool, version in software_versions.items():
        assert re.match(
            r"^\d.*|^[a-f0-9]{40}$", str(version)
        ), f"Version number for {tool} must start with a number, or be a Git SHA commit id. "
