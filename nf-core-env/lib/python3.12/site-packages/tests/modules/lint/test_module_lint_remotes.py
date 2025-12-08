import nf_core.modules.lint

from ...test_modules import TestModules
from ...utils import GITLAB_URL


class TestModulesLintRemotes(TestModules):
    """Test ModuleLint functionality with different remote sources"""

    def test_modules_lint_empty(self):
        """Test linting a pipeline with no modules installed"""
        self.mods_remove.remove("fastqc", force=True)
        self.mods_remove.remove("multiqc", force=True)
        nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        assert "No modules from https://github.com/nf-core/modules.git installed in pipeline" in self.caplog.text

    def test_modules_lint_no_gitlab(self):
        """Test linting a pipeline with no modules installed from gitlab"""
        self.mods_remove.remove("fastqc", force=True)
        self.mods_remove.remove("multiqc", force=True)
        nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir, remote_url=GITLAB_URL)
        assert f"No modules from {GITLAB_URL} installed in pipeline" in self.caplog.text

    def test_modules_lint_gitlab_modules(self):
        """Lint modules from a different remote"""
        self.mods_install_gitlab.install("fastqc")
        self.mods_install_gitlab.install("multiqc")
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir, remote_url=GITLAB_URL)
        module_lint.lint(print_results=False, all_modules=True)
        assert len(module_lint.failed) == 2
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0

    def test_modules_lint_multiple_remotes(self):
        """Lint modules from a different remote"""
        self.mods_install_gitlab.install("multiqc")
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir, remote_url=GITLAB_URL)
        module_lint.lint(print_results=False, all_modules=True)
        assert len(module_lint.failed) == 1
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0
