import json
from pathlib import Path

from git.repo import Repo

import nf_core.modules.lint
import nf_core.modules.patch

from ...test_modules import TestModules
from ...utils import GITLAB_NFTEST_BRANCH, GITLAB_URL


class TestModuleTests(TestModules):
    """Test module_tests.py functionality"""

    def test_modules_lint_snapshot_file(self):
        """Test linting a module with a snapshot file"""
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")
        assert len(module_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0

    def test_modules_lint_snapshot_file_missing_fail(self):
        """Test linting a module with a snapshot file missing, which should fail"""
        Path(
            self.nfcore_modules,
            "modules",
            "nf-core",
            "bpipe",
            "test",
            "tests",
            "main.nf.test.snap",
        ).unlink()
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")
        Path(
            self.nfcore_modules,
            "modules",
            "nf-core",
            "bpipe",
            "test",
            "tests",
            "main.nf.test.snap",
        ).touch()
        assert len(module_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0
        assert module_lint.failed[0].lint_test == "test_snapshot_exists"

    def test_modules_lint_snapshot_file_not_needed(self):
        """Test linting a module which doesn't need a snapshot file by removing the snapshot keyword in the main.nf.test file"""
        with open(
            Path(
                self.nfcore_modules,
                "modules",
                "nf-core",
                "bpipe",
                "test",
                "tests",
                "main.nf.test",
            )
        ) as fh:
            content = fh.read()
            new_content = content.replace("snapshot(", "snap (")
        with open(
            Path(
                self.nfcore_modules,
                "modules",
                "nf-core",
                "bpipe",
                "test",
                "tests",
                "main.nf.test",
            ),
            "w",
        ) as fh:
            fh.write(new_content)
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")
        assert len(module_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0

    def test_modules_missing_test_dir(self):
        """Test linting a module with a missing test directory"""
        Path(self.nfcore_modules, "modules", "nf-core", "bpipe", "test", "tests").rename(
            Path(self.nfcore_modules, "modules", "nf-core", "bpipe", "test", "tests.bak")
        )
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")
        Path(self.nfcore_modules, "modules", "nf-core", "bpipe", "test", "tests.bak").rename(
            Path(self.nfcore_modules, "modules", "nf-core", "bpipe", "test", "tests")
        )
        assert len(module_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) >= 0
        assert len(module_lint.warned) >= 0
        assert module_lint.failed[0].lint_test == "test_dir_exists"

    def test_modules_missing_test_main_nf(self):
        """Test linting a module with a missing test/main.nf file"""
        (self.bpipe_test_module_path / "tests" / "main.nf.test").rename(
            self.bpipe_test_module_path / "tests" / "main.nf.test.bak"
        )
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")
        (self.bpipe_test_module_path / "tests" / "main.nf.test.bak").rename(
            self.bpipe_test_module_path / "tests" / "main.nf.test"
        )
        assert len(module_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) >= 0
        assert len(module_lint.warned) >= 0
        assert module_lint.failed[0].lint_test == "test_main_nf_exists"

    def test_modules_unused_pytest_files(self):
        """Test linting a nf-test module with files still present in `tests/modules/`"""
        Path(self.nfcore_modules, "tests", "modules", "bpipe", "test").mkdir(parents=True, exist_ok=True)
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")
        Path(self.nfcore_modules, "tests", "modules", "bpipe", "test").rmdir()
        assert len(module_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) >= 0
        assert len(module_lint.warned) >= 0
        assert module_lint.failed[0].lint_test == "test_old_test_dir"

    def test_nftest_failing_linting(self):
        """Test linting a module which includes other modules in nf-test tests.
        Linting tests"""
        # Clone modules repo with testing modules
        tmp_dir = self.nfcore_modules.parent
        self.nfcore_modules = Path(tmp_dir, "modules-test")
        Repo.clone_from(GITLAB_URL, self.nfcore_modules, branch=GITLAB_NFTEST_BRANCH)

        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="kallisto/quant")

        assert len(module_lint.failed) == 2, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) >= 0
        assert len(module_lint.warned) >= 0
        assert module_lint.failed[0].lint_test == "meta_yml_valid"
        assert module_lint.failed[1].lint_test == "test_main_tags"
        assert "kallisto/index" in module_lint.failed[1].message

    def test_modules_absent_version(self):
        """Test linting a nf-test module if the versions is absent in the snapshot file `"""
        snap_file = self.bpipe_test_module_path / "tests" / "main.nf.test.snap"
        with open(snap_file) as fh:
            content = fh.read()
            new_content = content.replace("versions", "foo")
        with open(snap_file, "w") as fh:
            fh.write(new_content)
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")
        with open(snap_file, "w") as fh:
            fh.write(content)
        assert len(module_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) >= 0
        assert len(module_lint.warned) >= 0
        assert module_lint.failed[0].lint_test == "test_snap_versions"

    def test_modules_empty_file_in_snapshot(self):
        """Test linting a nf-test module with an empty file sha sum in the test snapshot, which should make it fail (if it is not a stub)"""
        snap_file = self.bpipe_test_module_path / "tests" / "main.nf.test.snap"
        snap = json.load(snap_file.open())
        content = snap_file.read_text()
        snap["my test"]["content"][0]["0"] = "test:md5,d41d8cd98f00b204e9800998ecf8427e"

        with open(snap_file, "w") as fh:
            json.dump(snap, fh)

        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")
        assert len(module_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0
        assert module_lint.failed[0].lint_test == "test_snap_md5sum"

        # reset the file
        with open(snap_file, "w") as fh:
            fh.write(content)

    def test_modules_empty_file_in_stub_snapshot(self):
        """Test linting a nf-test module with an empty file sha sum in the stub test snapshot, which should make it not fail"""
        snap_file = self.bpipe_test_module_path / "tests" / "main.nf.test.snap"
        snap = json.load(snap_file.open())
        content = snap_file.read_text()
        snap["my_test_stub"] = {"content": [{"0": "test:md5,d41d8cd98f00b204e9800998ecf8427e", "versions": {}}]}

        with open(snap_file, "w") as fh:
            json.dump(snap, fh)

        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")
        assert len(module_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0

        # Check for test_snap_md5sum in passed tests (handle both tuple and LintResult formats)
        found_test = False
        for x in module_lint.passed:
            test_name = None
            if hasattr(x, "lint_test"):
                test_name = x.lint_test
            elif isinstance(x, tuple) and len(x) > 0:
                test_name = x[0]

            if test_name == "test_snap_md5sum":
                found_test = True
                break

        assert found_test, "test_snap_md5sum not found in passed tests"

        # reset the file
        with open(snap_file, "w") as fh:
            fh.write(content)
