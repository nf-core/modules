import nf_core.modules.lint

from ...test_modules import TestModules


class TestModuleTodos(TestModules):
    """Test module_todos.py functionality"""

    def setUp(self):
        """Set up test fixtures by installing required modules"""
        super().setUp()
        # Install samtools/sort module for all tests in this class
        assert self.mods_install.install("samtools/sort")

    def test_module_todos_none(self):
        """Test module todos when no TODOs exist"""
        # Clean any TODO statements from files (they should be clean by default)
        module_dir = self.pipeline_dir / "modules" / "nf-core" / "samtools" / "sort"

        # Ensure main.nf has no TODO statements
        main_nf_path = module_dir / "main.nf"
        with open(main_nf_path) as fh:
            main_nf_content = fh.read()

        # Remove any TODO statements if they exist
        main_nf_content = main_nf_content.replace("TODO", "")
        with open(main_nf_path, "w") as fh:
            fh.write(main_nf_content)

        # Run lint on the module
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, module="samtools/sort", key=["module_todos"])

        # Should not have any warnings from TODOs
        warned_test_names = [test.lint_test for test in module_lint.warned]
        assert "module_todo" not in warned_test_names

    def test_module_todos_found_in_main_nf(self):
        """Test module todos when TODOs are found in main.nf"""
        # Add a TODO statement to main.nf
        module_dir = self.pipeline_dir / "modules" / "nf-core" / "samtools" / "sort"
        main_nf_path = module_dir / "main.nf"

        with open(main_nf_path, "a") as fh:
            fh.write("\n// TODO nf-core: This is a test TODO statement\n")

        # Run lint on the module
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, module="samtools/sort", key=["module_todos"])

        # Should have warning for TODO statement
        assert len(module_lint.warned) > 0, "Expected linting to warn due to TODO statement"
        warned_test_names = [test.lint_test for test in module_lint.warned]
        assert "module_todo" in warned_test_names

        # Check the specific warning message
        todo_warning = [test for test in module_lint.warned if test.lint_test == "module_todo"][0]
        assert "TODO" in todo_warning.message

    def test_module_todos_found_in_meta_yml(self):
        """Test module todos when TODOs are found in meta.yml"""
        # Add a TODO comment to meta.yml
        module_dir = self.pipeline_dir / "modules" / "nf-core" / "samtools" / "sort"
        meta_yml_path = module_dir / "meta.yml"

        with open(meta_yml_path, "a") as fh:
            fh.write("\n# TODO nf-core: Add more detailed description\n")

        # Run lint on the module
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, module="samtools/sort", key=["module_todos"])

        # Should have warning for TODO statement
        assert len(module_lint.warned) > 0, "Expected linting to warn due to TODO statement"
        warned_test_names = [test.lint_test for test in module_lint.warned]
        assert "module_todo" in warned_test_names

    def test_module_todos_multiple_found(self):
        """Test module todos when multiple TODOs are found"""
        # Add multiple TODO statements to different files
        module_dir = self.pipeline_dir / "modules" / "nf-core" / "samtools" / "sort"

        # Add TODO to main.nf
        main_nf_path = module_dir / "main.nf"
        with open(main_nf_path, "a") as fh:
            fh.write("\n// TODO nf-core: First TODO statement\n")
            fh.write("// TODO nf-core: Second TODO statement\n")

        # Add TODO to meta.yml
        meta_yml_path = module_dir / "meta.yml"
        with open(meta_yml_path, "a") as fh:
            fh.write("\n# TODO nf-core: Meta TODO statement\n")

        # Run lint on the module
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, module="samtools/sort", key=["module_todos"])

        # Should have multiple warnings for TODO statements
        todo_warnings = [test for test in module_lint.warned if test.lint_test == "module_todo"]
        assert len(todo_warnings) >= 3, f"Expected at least 3 TODO warnings, got {len(todo_warnings)}"
