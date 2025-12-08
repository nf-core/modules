import pytest

import nf_core.modules.lint

from ...test_modules import TestModules


class TestModuleChanges(TestModules):
    """Test module_changes.py functionality"""

    def setUp(self):
        """Set up test fixtures by installing required modules"""
        super().setUp()
        # Install samtools/sort module for all tests in this class
        assert self.mods_install.install("samtools/sort")

    def test_module_changes_unchanged(self):
        """Test module changes when module is unchanged"""
        # Run lint on the unchanged module
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, module="samtools/sort", key=["module_changes"])

        # Check that module_changes test passed (no changes detected)
        assert len(module_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"

        # Should have passed entries for files being up to date
        passed_test_names = [test.lint_test for test in module_lint.passed]
        assert "check_local_copy" in passed_test_names

    def test_module_changes_modified_main_nf(self):
        """Test module changes when main.nf is modified"""
        # Modify the main.nf file
        main_nf_path = self.pipeline_dir / "modules" / "nf-core" / "samtools" / "sort" / "main.nf"
        with open(main_nf_path, "a") as fh:
            fh.write("\n// This is a test modification\n")

        # Run lint on the modified module
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, module="samtools/sort", key=["module_changes"])

        # Check that module_changes test failed (changes detected)
        assert len(module_lint.failed) > 0, "Expected linting to fail due to modified file"

        # Should have failed entry for local copy not matching remote
        failed_test_names = [test.lint_test for test in module_lint.failed]
        assert "check_local_copy" in failed_test_names

    def test_module_changes_modified_meta_yml(self):
        """Test module changes when meta.yml is modified"""
        # Modify the meta.yml file
        meta_yml_path = self.pipeline_dir / "modules" / "nf-core" / "samtools" / "sort" / "meta.yml"
        with open(meta_yml_path, "a") as fh:
            fh.write("\n# This is a test comment\n")

        # Run lint on the modified module
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, module="samtools/sort", key=["module_changes"])

        # Check that module_changes test failed (changes detected)
        assert len(module_lint.failed) > 0, "Expected linting to fail due to modified file"

        # Should have failed entry for local copy not matching remote
        failed_test_names = [test.lint_test for test in module_lint.failed]
        assert "check_local_copy" in failed_test_names

    @pytest.mark.skip(reason="Patch testing requires complex setup - test framework needs improvement")
    def test_module_changes_patched_module(self):
        """Test module changes when module is patched"""
        # This test would require creating a patched module which is complex
        # in the current test framework. Skip for now until patch test infrastructure
        # is improved.
        pass
