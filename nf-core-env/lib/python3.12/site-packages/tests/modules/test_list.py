import json
from pathlib import Path

import yaml
from rich.console import Console

import nf_core.modules.list

from ..test_modules import TestModules
from ..utils import GITLAB_DEFAULT_BRANCH, GITLAB_URL


class TestModulesList(TestModules):
    def test_modules_list_remote(self):
        """Test listing available modules"""
        mods_list = nf_core.modules.list.ModuleList()
        listed_mods = mods_list.list_components()
        console = Console(record=True)
        console.print(listed_mods)
        output = console.export_text()
        assert "fastqc" in output

    def test_modules_list_remote_gitlab(self):
        """Test listing the modules in the remote gitlab repo"""
        mods_list = nf_core.modules.list.ModuleList(remote_url=GITLAB_URL, branch=GITLAB_DEFAULT_BRANCH)
        listed_mods = mods_list.list_components()
        console = Console(record=True)
        console.print(listed_mods)
        output = console.export_text()
        assert "fastqc" in output

    def test_modules_list_pipeline(self):
        """Test listing locally installed modules"""
        mods_list = nf_core.modules.list.ModuleList(self.pipeline_dir, remote=False)
        listed_mods = mods_list.list_components()
        console = Console(record=True)
        console.print(listed_mods)
        output = console.export_text()
        assert "fastqc" in output
        assert "multiqc" in output

    def test_modules_install_and_list_pipeline(self):
        """Test listing locally installed modules"""
        self.mods_install.install("trimgalore")
        mods_list = nf_core.modules.list.ModuleList(self.pipeline_dir, remote=False)
        listed_mods = mods_list.list_components()
        console = Console(record=True)
        console.print(listed_mods)
        output = console.export_text()
        assert "trimgalore" in output

    def test_modules_install_gitlab_and_list_pipeline(self):
        """Test listing locally installed modules"""
        self.mods_install_gitlab.install("fastqc")
        mods_list = nf_core.modules.list.ModuleList(self.pipeline_dir, remote=False)
        listed_mods = mods_list.list_components()
        console = Console(record=True)
        console.print(listed_mods)
        output = console.export_text()
        assert "fastqc" in output

    def test_modules_list_local_json(self):
        """Test listing locally installed modules as JSON"""
        mods_list = nf_core.modules.list.ModuleList(self.pipeline_dir, remote=False)
        listed_mods = str(mods_list.list_components(print_json=True))
        listed_mods = json.loads(listed_mods)
        assert "fastqc" in listed_mods
        assert "multiqc" in listed_mods

    def test_modules_list_remote_json(self) -> None:
        """Test listing available modules as JSON"""
        mods_list = nf_core.modules.list.ModuleList(remote=True)
        listed_mods: str = str(mods_list.list_components(print_json=True))
        listed_mods = json.loads(listed_mods)
        assert "fastqc" in listed_mods
        assert "multiqc" in listed_mods

    def test_modules_list_with_one_keyword(self):
        """Test listing available modules with one keyword"""
        mods_list = nf_core.modules.list.ModuleList(remote=True)
        listed_mods = mods_list.list_components(keywords=["qc"])
        console = Console(record=True)
        console.print(listed_mods)
        output = console.export_text()
        assert "multiqc" in output

    def test_modules_list_with_keywords(self):
        """Test listing available modules with multiple keywords"""
        mods_list = nf_core.modules.list.ModuleList(remote=True)
        listed_mods = mods_list.list_components(keywords=["fastq", "qc"])
        console = Console(record=True)
        console.print(listed_mods)
        output = console.export_text()
        assert "fastqc" in output

    def test_modules_list_with_unused_keyword(self):
        """Test listing available modules with an unused keyword"""
        mods_list = nf_core.modules.list.ModuleList(remote=True)
        with self.assertLogs(level="INFO") as log:
            listed_mods = mods_list.list_components(keywords=["you_will_never_find_me"])
            self.assertIn("No available", log.output[0])
            # expect empty list
            assert listed_mods == ""

    def test_modules_list_in_wrong_repo_fail(self):
        """Test listing available modules in a non-pipeline repo"""
        # modify repotype in .nf-core.yml
        with open(Path(self.pipeline_dir, ".nf-core.yml")) as fh:
            nf_core_yml = yaml.safe_load(fh)
        nf_core_yml_orig = nf_core_yml.copy()
        nf_core_yml["repository_type"] = "modules"
        nf_core_yml["org_path"] = "nf-core"

        print(nf_core_yml)
        with open(Path(self.pipeline_dir, ".nf-core.yml"), "w") as fh:
            yaml.safe_dump(nf_core_yml, fh)
        # expect error logged
        with self.assertLogs(level="ERROR") as log:
            mods_list = nf_core.modules.list.ModuleList(self.pipeline_dir, remote=False)
            listed_mods = mods_list.list_components()
            self.assertIn("must be run from a pipeline directory", log.output[0])
            # expect empty list
            assert listed_mods == ""
        # restore .nf-core.yml
        with open(Path(self.pipeline_dir, ".nf-core.yml"), "w") as fh:
            yaml.safe_dump(nf_core_yml_orig, fh)
