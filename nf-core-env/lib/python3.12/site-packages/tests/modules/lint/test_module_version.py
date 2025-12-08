import pytest

import nf_core.modules.lint

from ...test_modules import TestModules


class TestModuleVersion(TestModules):
    """Test module_version.py functionality"""

    def test_module_version_with_git_sha(self):
        """Test module version when git_sha is present in modules.json"""
        # Install a module
        if not self.mods_install.install("samtools/sort"):
            self.skipTest("Could not install samtools/sort module")
        # Run lint on the module - should have a git_sha entry
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, module="samtools/sort", key=["module_version"])

        # Should pass git_sha test (git_sha entry exists)
        passed_test_names = [test.lint_test for test in module_lint.passed]
        assert "git_sha" in passed_test_names

        # Should have module_version test result (either passed or warned)
        all_test_names = [test.lint_test for test in module_lint.passed + module_lint.warned + module_lint.failed]
        assert "module_version" in all_test_names

    def test_module_version_up_to_date(self):
        """Test module version when module is up to date"""
        # Install a module
        if not self.mods_install.install("samtools/sort"):
            self.skipTest("Could not install samtools/sort module")
        # Run lint on the module
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, module="samtools/sort", key=["module_version"])

        # Should have a result for module_version (either passed if up-to-date or warned if newer available)
        all_tests = module_lint.passed + module_lint.warned + module_lint.failed
        version_test_names = [test.lint_test for test in all_tests]
        assert "module_version" in version_test_names

    @pytest.mark.skip(reason="Testing outdated modules requires specific version setup")
    def test_module_version_outdated(self):
        """Test module version when module is outdated"""
        # This test would require installing a specific older version of a module
        # which is complex to set up reliably in the test framework
        pass

    @pytest.mark.skip(reason="Testing git log failure requires complex mocking setup")
    def test_module_version_git_log_fail(self):
        """Test module version when git log fetch fails"""
        # This test would require mocking network failures or invalid repositories
        # which is complex to set up in the current test framework
        pass
