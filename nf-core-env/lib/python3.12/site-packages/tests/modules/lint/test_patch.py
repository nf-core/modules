from pathlib import Path

import nf_core.modules.install
import nf_core.modules.lint
import nf_core.modules.patch
from nf_core.utils import set_wd

from ...test_modules import TestModules
from ...utils import GITLAB_URL
from ..test_patch import BISMARK_ALIGN, CORRECT_SHA, PATCH_BRANCH, REPO_NAME, modify_main_nf


class TestPatch(TestModules):
    """Test patch.py functionality"""

    def _setup_patch(self, pipeline_dir: str | Path, modify_module: bool):
        install_obj = nf_core.modules.install.ModuleInstall(
            pipeline_dir,
            prompt=False,
            force=False,
            remote_url=GITLAB_URL,
            branch=PATCH_BRANCH,
            sha=CORRECT_SHA,
        )

        # Install the module
        install_obj.install(BISMARK_ALIGN)

        if modify_module:
            # Modify the module
            module_path = Path(pipeline_dir, "modules", REPO_NAME, BISMARK_ALIGN)
            modify_main_nf(module_path / "main.nf")

    def test_modules_lint_patched_modules(self):
        """
        Test creating a patch file and applying it to a new version of the the files
        """
        self._setup_patch(str(self.pipeline_dir), True)

        # Create a patch file
        patch_obj = nf_core.modules.patch.ModulePatch(self.pipeline_dir, GITLAB_URL, PATCH_BRANCH)
        patch_obj.patch(BISMARK_ALIGN)

        # change temporarily working directory to the pipeline directory
        # to avoid error from try_apply_patch() during linting
        with set_wd(self.pipeline_dir):
            module_lint = nf_core.modules.lint.ModuleLint(
                directory=self.pipeline_dir,
                remote_url=GITLAB_URL,
                branch=PATCH_BRANCH,
                hide_progress=True,
            )
            module_lint.lint(
                print_results=False,
                all_modules=True,
            )

        assert len(module_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0
