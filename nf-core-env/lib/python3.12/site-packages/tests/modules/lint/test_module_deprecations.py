import nf_core.modules.lint

from ...test_modules import TestModules


class TestModuleDeprecations(TestModules):
    """Test module_deprecations.py functionality"""

    def setUp(self):
        """Set up test fixtures by installing required modules"""
        super().setUp()
        # Install samtools/sort module for all tests in this class
        assert self.mods_install.install("samtools/sort")

    def test_module_deprecations_none(self):
        """Test module deprecations when no deprecations exist"""
        # Run lint on the module
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, module="samtools/sort", key=["module_deprecations"])

        # Should not have any failures from deprecations
        failed_test_names = [test.lint_test for test in module_lint.failed]
        assert "module_deprecations" not in failed_test_names

    def test_module_deprecations_functions_nf(self):
        """Test module deprecations when functions.nf exists"""
        # Create a deprecated functions.nf file
        module_dir = self.pipeline_dir / "modules" / "nf-core" / "samtools" / "sort"
        functions_nf_path = module_dir / "functions.nf"

        # Create the deprecated functions.nf file
        with open(functions_nf_path, "w") as fh:
            fh.write("// Deprecated functions.nf file\n")

        # Run lint on the module
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, module="samtools/sort", key=["module_deprecations"])

        # Should have failure for deprecated functions.nf file
        assert len(module_lint.failed) > 0, "Expected linting to fail due to deprecated functions.nf file"
        failed_test_names = [test.lint_test for test in module_lint.failed]
        assert "module_deprecations" in failed_test_names

        # Check the specific failure message
        deprecation_failure = [test for test in module_lint.failed if test.lint_test == "module_deprecations"][0]
        assert "functions.nf" in deprecation_failure.message
        assert "Deprecated" in deprecation_failure.message

    def test_module_deprecations_no_functions_nf(self):
        """Test module deprecations when no functions.nf exists"""
        # Ensure no functions.nf file exists (should be default)
        module_dir = self.pipeline_dir / "modules" / "nf-core" / "samtools" / "sort"
        functions_nf_path = module_dir / "functions.nf"

        # Remove functions.nf if it somehow exists
        if functions_nf_path.exists():
            functions_nf_path.unlink()

        # Run lint on the module
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, module="samtools/sort", key=["module_deprecations"])

        # Should not have any failures from deprecations
        failed_test_names = [test.lint_test for test in module_lint.failed]
        assert "module_deprecations" not in failed_test_names
