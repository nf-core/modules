from pathlib import Path

import yaml


def actions_awsfulltest(self) -> dict[str, list[str]]:
    """Checks the GitHub Actions awsfulltest is valid.

    In addition to small test datasets run on GitHub Actions, we provide the possibility of testing the pipeline on full size datasets on AWS.
    This should ensure that the pipeline runs as expected on AWS and provide a resource estimation.

    The GitHub Actions workflow is called ``awsfulltest.yml``, and it can be found in the ``.github/workflows/`` directory.

    .. warning::  This workflow incurs AWS costs, therefore it should only be triggered for pipeline releases:
                  ``release`` (after the pipeline release) and ``workflow_dispatch``.

    .. note::  You can manually trigger the AWS tests by going to the `Actions` tab on the pipeline GitHub repository and selecting the
                  `nf-core AWS full size tests` workflow on the left.

    .. tip::  For tests on full data prior to release, `Seqera Platform <https://seqera.io/platform/>`_ launch feature can be employed.

    The ``.github/workflows/awsfulltest.yml`` file is tested for the following:

    * Must be turned on ``workflow_dispatch``.
    * Must be turned on for ``release`` with ``types: [published]``
    * Should run the profile ``test_full`` that should be edited to provide the links to full-size datasets. If it runs the profile ``test``, a warning is given.
    """
    passed = []
    warned = []
    failed = []

    fn = Path(self.wf_path, ".github", "workflows", "awsfulltest.yml")
    if fn.is_file():
        try:
            with open(fn) as fh:
                wf = yaml.safe_load(fh)
        except Exception as e:
            return {"failed": [f"Could not parse yaml file: {fn}, {e}"]}

        aws_profile = "-profile test "

        # Check that the action is only turned on for published releases
        try:
            if wf[True]["pull_request_review"]["types"] != ["submitted"]:
                raise AssertionError()
            if "workflow_dispatch" not in wf[True]:
                raise AssertionError()
            if wf[True]["release"]["types"] != ["published"]:
                raise AssertionError()
        except (AssertionError, KeyError, TypeError):
            failed.append("`.github/workflows/awsfulltest.yml` is not triggered correctly")
        else:
            passed.append("`.github/workflows/awsfulltest.yml` is triggered correctly")

        # Warn if `-profile test` is still unchanged
        try:
            steps = wf["jobs"]["run-platform"]["steps"]
            if not any(aws_profile in step["run"] for step in steps if "run" in step.keys()):
                raise AssertionError()
        except (AssertionError, KeyError, TypeError):
            passed.append("`.github/workflows/awsfulltest.yml` does not use `-profile test`")
        else:
            warned.append("`.github/workflows/awsfulltest.yml` should test full datasets, not `-profile test`")

    return {"passed": passed, "warned": warned, "failed": failed}
