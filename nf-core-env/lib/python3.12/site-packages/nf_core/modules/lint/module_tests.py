"""
Lint the tests of a module in nf-core/modules
"""

import json
import logging
import re
from pathlib import Path

from nf_core.components.lint import LintExceptionError
from nf_core.components.nfcore_component import NFCoreComponent

log = logging.getLogger(__name__)


def module_tests(_, module: NFCoreComponent, allow_missing: bool = False):
    """
    Lint the tests of a module in ``nf-core/modules``

    It verifies that the test directory exists
    and contains a ``main.nf.test`` and a ``main.nf.test.snap``

    """
    if module.nftest_testdir is None:
        if allow_missing:
            module.warned.append(
                (
                    "module_tests",
                    "test_dir_exists",
                    "nf-test directory is missing",
                    Path(module.component_dir, "tests"),
                )
            )
            return
        raise LintExceptionError("Module does not have a `tests` dir")

    if module.nftest_main_nf is None:
        if allow_missing:
            module.warned.append(
                (
                    "module_tests",
                    "test_main_nf_exists",
                    "test `main.nf.test` does not exist",
                    Path(module.component_dir, "tests", "main.nf.test"),
                )
            )
            return
        raise LintExceptionError("Module does not have a `tests` dir")

    repo_dir = module.component_dir.parts[: module.component_dir.parts.index(module.component_name.split("/")[0])][-1]
    test_dir = Path(module.base_dir, "tests", "modules", repo_dir, module.component_name)
    pytest_main_nf = Path(test_dir, "main.nf")
    is_pytest = pytest_main_nf.is_file()
    if module.nftest_testdir.is_dir():
        module.passed.append(
            ("module_tests", "test_dir_exists", "nf-test test directory exists", module.nftest_testdir)
        )
    else:
        if is_pytest:
            module.warned.append(
                (
                    "module_tests",
                    "test_dir_exists",
                    "nf-test directory is missing",
                    module.nftest_testdir,
                )
            )
        else:
            module.failed.append(
                (
                    "module_tests",
                    "test_dir_exists",
                    "nf-test directory is missing",
                    module.nftest_testdir,
                )
            )
        return

    # Lint the test main.nf file
    if module.nftest_main_nf.is_file():
        module.passed.append(
            ("module_tests", "test_main_nf_exists", "test `main.nf.test` exists", module.nftest_main_nf)
        )
    else:
        if is_pytest:
            module.warned.append(
                (
                    "module_tests",
                    "test_main_nf_exists",
                    "test `main.nf.test` does not exist",
                    module.nftest_main_nf,
                )
            )
        else:
            module.failed.append(
                (
                    "module_tests",
                    "test_main_nf_exists",
                    "test `main.nf.test` does not exist",
                    module.nftest_main_nf,
                )
            )

    if module.nftest_main_nf.is_file():
        # Check if main.nf.test.snap file exists, if 'snap(' is inside main.nf.test
        with open(module.nftest_main_nf) as fh:
            if "snapshot(" in fh.read():
                snap_file = module.nftest_testdir / "main.nf.test.snap"

                if snap_file.is_file():
                    module.passed.append(
                        (
                            "module_tests",
                            "test_snapshot_exists",
                            "snapshot file `main.nf.test.snap` exists",
                            snap_file,
                        )
                    )
                    # Validate no empty files
                    with open(snap_file) as snap_fh:
                        try:
                            snap_content = json.load(snap_fh)
                            for test_name in snap_content.keys():
                                if "d41d8cd98f00b204e9800998ecf8427e" in str(snap_content[test_name]):
                                    if "stub" not in test_name:
                                        module.failed.append(
                                            (
                                                "module_tests",
                                                "test_snap_md5sum",
                                                "md5sum for empty file found: d41d8cd98f00b204e9800998ecf8427e",
                                                snap_file,
                                            )
                                        )
                                    else:
                                        module.passed.append(
                                            (
                                                "module_tests",
                                                "test_snap_md5sum",
                                                "md5sum for empty file found, but it is a stub test",
                                                snap_file,
                                            )
                                        )
                                else:
                                    module.passed.append(
                                        (
                                            "module_tests",
                                            "test_snap_md5sum",
                                            "no md5sum for empty file found",
                                            snap_file,
                                        )
                                    )
                                if "7029066c27ac6f5ef18d660d5741979a" in str(snap_content[test_name]):
                                    if "stub" not in test_name:
                                        module.failed.append(
                                            (
                                                "module_tests",
                                                "test_snap_md5sum",
                                                "md5sum for compressed empty file found: 7029066c27ac6f5ef18d660d5741979a",
                                                snap_file,
                                            )
                                        )
                                    else:
                                        module.passed.append(
                                            (
                                                "module_tests",
                                                "test_snap_md5sum",
                                                "md5sum for compressed empty file found, but it is a stub test",
                                                snap_file,
                                            )
                                        )
                                else:
                                    module.passed.append(
                                        (
                                            "module_tests",
                                            "test_snap_md5sum",
                                            "no md5sum for compressed empty file found",
                                            snap_file,
                                        )
                                    )
                            if "versions" in str(snap_content[test_name]) or "versions" in str(snap_content.keys()):
                                module.passed.append(
                                    (
                                        "module_tests",
                                        "test_snap_versions",
                                        "versions found in snapshot file",
                                        snap_file,
                                    )
                                )
                            else:
                                module.failed.append(
                                    (
                                        "module_tests",
                                        "test_snap_versions",
                                        "versions not found in snapshot file",
                                        snap_file,
                                    )
                                )
                        except json.decoder.JSONDecodeError as e:
                            module.failed.append(
                                (
                                    "module_tests",
                                    "test_snapshot_exists",
                                    f"snapshot file `main.nf.test.snap` can't be read: {e}",
                                    snap_file,
                                )
                            )
                else:
                    module.failed.append(
                        (
                            "module_tests",
                            "test_snapshot_exists",
                            "snapshot file `main.nf.test.snap` does not exist",
                            snap_file,
                        )
                    )
            # Verify that tags are correct.
            main_nf_tags = module._get_main_nf_tags(module.nftest_main_nf)
            not_alphabet = re.compile(r"[^a-zA-Z]")
            org_alp = not_alphabet.sub("", module.org)
            org_alphabet = org_alp if org_alp != "" else "nfcore"
            required_tags = ["modules", f"modules_{org_alphabet}", module.component_name]
            if module.component_name.count("/") == 1:
                required_tags.append(module.component_name.split("/")[0])
            chained_components_tags = module._get_included_components_in_chained_tests(module.nftest_main_nf)
            missing_tags = []
            log.debug(f"Required tags: {required_tags}")
            log.debug(f"Included components for chained nf-tests: {chained_components_tags}")
            for tag in set(required_tags + chained_components_tags):
                if tag not in main_nf_tags:
                    missing_tags.append(tag)
            if len(missing_tags) == 0:
                module.passed.append(
                    (
                        "module_tests",
                        "test_main_tags",
                        "Tags adhere to guidelines",
                        module.nftest_main_nf,
                    )
                )
            else:
                module.failed.append(
                    (
                        "module_tests",
                        "test_main_tags",
                        f"Tags do not adhere to guidelines. Tags missing in `main.nf.test`: `{','.join(missing_tags)}`",
                        module.nftest_main_nf,
                    )
                )

    # Check that the old test directory does not exist
    if not is_pytest:
        old_test_dir = Path(module.base_dir, "tests", "modules", module.component_name)
        if old_test_dir.is_dir():
            module.failed.append(
                (
                    "module_tests",
                    "test_old_test_dir",
                    f"Pytest files are still present at `{Path('tests', 'modules', module.component_name)}`. Please remove this directory and its contents.",
                    old_test_dir,
                )
            )
        else:
            module.passed.append(
                (
                    "module_tests",
                    "test_old_test_dir",
                    "Old pytests don't exist for this module",
                    old_test_dir,
                )
            )
