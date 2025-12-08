import logging
import re
from pathlib import Path

from nf_core.utils import load_tools_config

log = logging.getLogger(__name__)


def nf_test_content(self) -> dict[str, list[str]]:
    """Checks that the pipeline nf-test files have the appropriate content.

    This lint test checks the following files and content of these files:

    * ``*.nf.test`` files should specify the ``outdir`` parameter:

    .. code-block:: groovy

        when {
            params {
                outdir = "$outputDir"
            }
        }

    * A `versions.yml` file should be included in the snapshot of all `*.nf.test` files

    * The `nextflow.config` file should contain:

        .. code-block:: groovy

            modules_testdata_base_path = <path>

        and

        .. code-block:: groovy

            pipelines_testdata_base_path = <path>

        And should set the correct resource limits, as defined in the `test` profile

    * The `nf-test.config` file should:
        * Make sure tests are relative to root directory

        .. code-block:: groovy

            testsDir "."

        * Ensure a user-configurable nf-test directory

        .. code-block:: groovy

            workDir System.getenv("NFT_WORKDIR") ?: ".nf-test"

        * Use a test specific config

        .. code-block:: groovy

            configFile "tests/nextflow.config"

    All these checks can be skipped in the `.nf-core.yml` file using:

    .. code-block:: yaml

        lint:
            nf_test_content: False

    or

    .. code-block:: yaml

        lint:
            nf_test_content:
                - tests/<test_name>.nf.test
                - tests/nextflow.config
                - nf-test.config
    """
    passed: list[str] = []
    failed: list[str] = []
    ignored: list[str] = []

    _, pipeline_conf = load_tools_config(self.wf_path)
    lint_conf = getattr(pipeline_conf, "lint", None) or None
    nf_test_content_conf = getattr(lint_conf, "nf_test_content", None) or None

    # Content of *.nf.test files
    test_fns = list(Path(self.wf_path, "tests").glob("*.nf.test"))
    test_checks: dict[str, dict[str, str | bool]] = {
        "outdir": {
            "pattern": r"outdir *= *[\"']\${?outputDir}?[\"']",
            "description": "`outdir` parameter",
            "failure_msg": 'does not contain `outdir` parameter, it should contain `outdir = "$outputDir"`',
            "when_block": True,
        },
        "versions.yml": {
            "pattern": r"versions\.yml",
            "description": "snapshots a 'versions.yml' file",
            "failure_msg": "does not snapshot a 'versions.yml' file",
            "when_block": False,
        },
    }

    for test_fn in test_fns:
        if nf_test_content_conf is not None and (
            not nf_test_content_conf or str(test_fn.relative_to(self.wf_path)) in nf_test_content_conf
        ):
            ignored.append(f"'{test_fn.relative_to(self.wf_path)}' checking ignored")
            continue

        checks_passed = {check: False for check in test_checks}
        with open(test_fn) as fh:
            for line in fh:
                for check_name, check_info in test_checks.items():
                    if check_info["when_block"] and "when" in line:
                        while "}\n" not in line:
                            line = next(fh)
                            if re.search(str(check_info["pattern"]), line):
                                passed.append(
                                    f"'{test_fn.relative_to(self.wf_path)}' contains {check_info['description']}"
                                )
                                checks_passed[check_name] = True
                                break
                    elif not check_info["when_block"] and re.search(str(check_info["pattern"]), line):
                        passed.append(f"'{test_fn.relative_to(self.wf_path)}' {check_info['description']}")
                        checks_passed[check_name] = True

        for check_name, check_info in test_checks.items():
            if not checks_passed[check_name]:
                failed.append(f"'{test_fn.relative_to(self.wf_path)}' {check_info['failure_msg']}")

    # Content of nextflow.config file
    conf_fn = Path(self.wf_path, "tests", "nextflow.config")

    config_checks: dict[str, dict[str, str]] = {
        "modules_testdata_base_path": {
            "pattern": "modules_testdata_base_path",
            "description": "`modules_testdata_base_path`",
        },
        "pipelines_testdata_base_path": {
            "pattern": "pipelines_testdata_base_path",
            "description": "`pipelines_testdata_base_path`",
        },
    }
    # Check if tests/nextflow.config is present
    if not conf_fn.exists():
        failed.append(f"'{conf_fn.relative_to(self.wf_path)}' does not exist")
    else:
        if nf_test_content_conf is None or str(conf_fn.relative_to(self.wf_path)) not in nf_test_content_conf:
            checks_passed = {check: False for check in config_checks}
            with open(conf_fn) as fh:
                for line in fh:
                    line = line.strip()
                    for check_name, config_check_info in config_checks.items():
                        if re.search(str(config_check_info["pattern"]), line):
                            passed.append(
                                f"'{conf_fn.relative_to(self.wf_path)}' contains {config_check_info['description']}"
                            )
                            checks_passed[check_name] = True
            for check_name, config_check_info in config_checks.items():
                if not checks_passed[check_name]:
                    failed.append(
                        f"'{conf_fn.relative_to(self.wf_path)}' does not contain {config_check_info['description']}"
                    )
        else:
            ignored.append(f"'{conf_fn.relative_to(self.wf_path)}' checking ignored")

    # Content of nf-test.config file
    nf_test_conf_fn = Path(self.wf_path, "nf-test.config")
    nf_test_checks: dict[str, dict[str, str]] = {
        "testsDir": {
            "pattern": r'testsDir "\."',
            "description": "sets a `testsDir`",
            "failure_msg": 'does not set a `testsDir`, it should contain `testsDir "."`',
        },
        "workDir": {
            "pattern": r'workDir System\.getenv\("NFT_WORKDIR"\) \?: "\.nf-test"',
            "description": "sets a `workDir`",
            "failure_msg": 'does not set a `workDir`, it should contain `workDir System.getenv("NFT_WORKDIR") ?: ".nf-test"`',
        },
        "configFile": {
            "pattern": r'configFile "tests/nextflow\.config"',
            "description": "sets a `configFile`",
            "failure_msg": 'does not set a `configFile`, it should contain `configFile "tests/nextflow.config"`',
        },
    }

    if not nf_test_conf_fn.exists():
        failed.append(f"'{nf_test_conf_fn.relative_to(self.wf_path)}' does not exist")
    else:
        if nf_test_content_conf is None or str(nf_test_conf_fn.relative_to(self.wf_path)) not in nf_test_content_conf:
            checks_passed = {check: False for check in nf_test_checks}
            with open(nf_test_conf_fn) as fh:
                for line in fh:
                    line = line.strip()
                    for check_name, nf_test_check_info in nf_test_checks.items():
                        if re.search(str(nf_test_check_info["pattern"]), line):
                            passed.append(
                                f"'{nf_test_conf_fn.relative_to(self.wf_path)}' {nf_test_check_info['description']}"
                            )
                            checks_passed[check_name] = True
            for check_name, nf_test_check_info in nf_test_checks.items():
                if not checks_passed[check_name]:
                    failed.append(f"'{nf_test_conf_fn.relative_to(self.wf_path)}' {nf_test_check_info['failure_msg']}")
        else:
            ignored.append(f"'{nf_test_conf_fn.relative_to(self.wf_path)}' checking ignored")

    return {"passed": passed, "failed": failed, "ignored": ignored}
