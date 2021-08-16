from pathlib import Path
import pytest
import yaml


def _get_workflow_names():
    """Get all names of all workflows which have a test.yml in the tests directory.

    To do so, recursively finds all test.yml files and parses their content.
    """
    here = Path(__file__).parent.resolve()
    pytest_workflow_files = here.glob("**/test.yml")
    for f in pytest_workflow_files:
        test_config = yaml.safe_load(f.read_text())
        for workflow in test_config:
            yield workflow["name"]


@pytest.mark.workflow(*_get_workflow_names())
def test_ensure_valid_version_yml(workflow_dir):
    workflow_dir = Path(workflow_dir)
    software_name = workflow_dir.name.split("_")[0].lower()
    versions_yml = (workflow_dir / f"output/{software_name}/versions.yml").read_text()
    assert "END_VERSIONS" not in versions_yml
    # Raises an exception if yaml is not valid
    yaml.safe_load(versions_yml)
