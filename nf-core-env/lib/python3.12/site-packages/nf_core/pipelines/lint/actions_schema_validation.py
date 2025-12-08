import logging
from pathlib import Path
from typing import Any

import jsonschema
import requests
import yaml


def actions_schema_validation(self) -> dict[str, list[str]]:
    """Checks that the GitHub Action workflow yml/yaml files adhere to the correct schema

    nf-core pipelines use GitHub actions workflows to run CI tests, check formatting and also linting, among others.
    These workflows are defined by ``yml`` scripts in ``.github/workflows/``. This lint test verifies that these scripts are valid
    by comparing them against the `JSON schema for GitHub workflows <https://www.schemastore.org/github-workflow>`_.

    To pass this test, make sure that all your workflows contain the required properties ``on`` and ``jobs`` and that
    all other properties are of the correct type, as specified in the schema (link above).
    """
    passed: list[str] = []
    failed: list[str] = []
    warned: list[str] = []

    # Only show error messages from schema
    logging.getLogger("nf_core.pipelines.schema").setLevel(logging.ERROR)

    # Get all workflow files
    action_workflows = list(Path(self.wf_path).glob(".github/workflows/*.y*ml"))

    # Load the GitHub workflow schema
    r = requests.get("https://www.schemastore.org/github-workflow", allow_redirects=True)
    # handle "Service Unavailable" error
    if r.status_code not in [200, 301]:
        warned.append(
            f"Failed to fetch schema: Response code for `https://www.schemastore.org/github-workflow` was {r.status_code}"
        )
        return {"passed": passed, "failed": failed, "warned": warned}
    schema: dict[str, Any] = r.json()

    # Validate all workflows against the schema
    for wf_path in action_workflows:
        wf = wf_path.name

        # load workflow
        try:
            with open(wf_path) as fh:
                wf_json = yaml.safe_load(fh)
        except Exception as e:
            failed.append(f"Could not parse yaml file: {wf}, {e}")
            continue

        # yaml parses 'on' as True --> try to fix it before schema validation
        try:
            wf_json["on"] = wf_json.pop(True)
        except Exception:
            failed.append(f"Missing 'on' keyword in {wf}")

        # Validate the workflow
        try:
            jsonschema.validate(wf_json, schema)
            passed.append(f"Workflow validation passed: {wf}")
        except Exception as e:
            failed.append(f"Workflow validation failed for {wf}: {e}")

    return {"passed": passed, "failed": failed, "warned": warned}
