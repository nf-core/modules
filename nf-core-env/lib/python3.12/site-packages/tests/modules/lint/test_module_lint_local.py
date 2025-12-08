import shutil
from pathlib import Path

import nf_core.modules.lint

from ...test_modules import TestModules


class TestModulesLintLocal(TestModules):
    """Test ModuleLint functionality with local modules"""

    def setUp(self):
        """Set up test fixtures by installing required modules"""
        super().setUp()
        # Install trimgalore module for all tests in this class
        if not self.mods_install.install("trimgalore"):
            self.skipTest("Could not install trimgalore module")

    def test_modules_lint_local(self):
        """Test linting local modules"""
        installed = Path(self.pipeline_dir, "modules", "nf-core", "trimgalore")
        local = Path(self.pipeline_dir, "modules", "local", "trimgalore")
        shutil.move(installed, local)
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, local=True)
        assert len(module_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0

    def test_modules_lint_local_missing_files(self):
        """Test linting local modules with missing files"""
        installed = Path(self.pipeline_dir, "modules", "nf-core", "trimgalore")
        local = Path(self.pipeline_dir, "modules", "local", "trimgalore")
        shutil.move(installed, local)
        Path(self.pipeline_dir, "modules", "local", "trimgalore", "environment.yml").unlink()
        Path(self.pipeline_dir, "modules", "local", "trimgalore", "meta.yml").unlink()
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, local=True)
        assert len(module_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0
        warnings = [x.message for x in module_lint.warned]
        assert "Module's `environment.yml` does not exist" in warnings
        assert "Module `meta.yml` does not exist" in warnings

    def test_modules_lint_local_old_format(self):
        """Test linting local modules in old format"""
        Path(self.pipeline_dir, "modules", "local").mkdir()
        installed = Path(self.pipeline_dir, "modules", "nf-core", "trimgalore", "main.nf")
        local = Path(self.pipeline_dir, "modules", "local", "trimgalore.nf")
        shutil.move(installed, local)
        self.mods_remove.remove("trimgalore", force=True)
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, local=True)
        assert len(module_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0
