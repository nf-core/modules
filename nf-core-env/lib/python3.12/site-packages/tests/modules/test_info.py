from rich.console import Console

import nf_core.modules.info

from ..test_modules import TestModules
from ..utils import GITLAB_DEFAULT_BRANCH, GITLAB_URL


class TestModulesCreate(TestModules):
    def test_modules_info_remote(self):
        """Test getting info about a remote module"""
        mods_info = nf_core.modules.info.ModuleInfo(self.pipeline_dir, "fastqc")
        mods_info_output = mods_info.get_component_info()
        console = Console(record=True)
        console.print(mods_info_output)
        output = console.export_text()

        assert "Module: fastqc" in output
        assert "Inputs" in output
        assert "Outputs" in output

    def test_modules_info_remote_gitlab(self):
        """Test getting info about a module in the remote gitlab repo"""
        mods_info = nf_core.modules.info.ModuleInfo(
            self.pipeline_dir, "fastqc", remote_url=GITLAB_URL, branch=GITLAB_DEFAULT_BRANCH
        )
        mods_info_output = mods_info.get_component_info()
        console = Console(record=True)
        console.print(mods_info_output)
        output = console.export_text()

        assert "Module: fastqc" in output
        assert "Inputs" in output
        assert "Outputs" in output
        assert "--git-remote" in output

    def test_modules_info_local(self):
        """Test getting info about a locally installed module"""
        self.mods_install.install("trimgalore")
        mods_info = nf_core.modules.info.ModuleInfo(self.pipeline_dir, "trimgalore")
        mods_info_output = mods_info.get_component_info()
        console = Console(record=True)
        console.print(mods_info_output)
        output = console.export_text()

        assert "Module: trimgalore" in output
        assert "Inputs" in output
        assert "Outputs" in output
        assert "Location" in output

    def test_modules_info_in_modules_repo(self):
        """Test getting info about a module in the modules repo"""
        mods_info = nf_core.modules.info.ModuleInfo(self.nfcore_modules, "fastqc")
        mods_info.local = True
        mods_info_output = mods_info.get_component_info()
        console = Console(record=True)
        console.print(mods_info_output)
        output = console.export_text()

        assert "Module: fastqc" in output
        assert "Inputs" in output
        assert "Outputs" in output
