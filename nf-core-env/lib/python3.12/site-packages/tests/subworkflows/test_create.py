import shutil
from pathlib import Path
from unittest import mock

import pytest
import yaml
from git.repo import Repo

import nf_core.subworkflows
from tests.utils import GITLAB_SUBWORKFLOWS_ORG_PATH_BRANCH, GITLAB_URL

from ..test_subworkflows import TestSubworkflows


class TestSubworkflowsCreate(TestSubworkflows):
    def test_subworkflows_create_succeed(self):
        """Succeed at creating a subworkflow from the template inside a pipeline"""
        subworkflow_create = nf_core.subworkflows.SubworkflowCreate(
            self.pipeline_dir, "test_subworkflow_local", "@author", True
        )
        subworkflow_create.create()
        assert Path(self.pipeline_dir, "subworkflows", "local", "test_subworkflow_local/main.nf").exists()

    def test_subworkflows_create_fail_exists(self):
        """Fail at creating the same subworkflow twice"""
        subworkflow_create = nf_core.subworkflows.SubworkflowCreate(
            self.pipeline_dir, "test_subworkflow2", "@author", False
        )
        subworkflow_create.create()
        with pytest.raises(UserWarning) as excinfo:
            subworkflow_create.create()
        assert "subworkflow directory exists" in str(excinfo.value)

    def test_subworkflows_create_nfcore_modules(self):
        """Create a subworkflow in nf-core/modules clone"""
        subworkflow_create = nf_core.subworkflows.SubworkflowCreate(
            self.nfcore_modules, "test_subworkflow", "@author", force=True
        )
        subworkflow_create.create()
        assert Path(self.nfcore_modules, "subworkflows", "nf-core", "test_subworkflow", "main.nf").exists()

        assert Path(
            self.nfcore_modules, "subworkflows", "nf-core", "test_subworkflow", "tests", "main.nf.test"
        ).exists()

    @mock.patch("rich.prompt.Confirm.ask")
    def test_subworkflows_migrate(self, mock_rich_ask):
        """Create a subworkflow with the --migrate-pytest option to convert pytest to nf-test"""
        pytest_dir = Path(self.nfcore_modules, "tests", "subworkflows", "nf-core", "bam_stats_samtools")
        subworkflow_dir = Path(self.nfcore_modules, "subworkflows", "nf-core", "bam_stats_samtools")

        # Clone modules repo with pytests
        shutil.rmtree(self.nfcore_modules)
        Repo.clone_from(GITLAB_URL, self.nfcore_modules, branch=GITLAB_SUBWORKFLOWS_ORG_PATH_BRANCH)
        with open(subworkflow_dir / "main.nf") as fh:
            old_main_nf = fh.read()
        with open(subworkflow_dir / "meta.yml") as fh:
            old_meta_yml = fh.read()

        # Create a subworkflow with --migrate-pytest
        mock_rich_ask.return_value = True
        subworkflow_create = nf_core.subworkflows.SubworkflowCreate(
            self.nfcore_modules, "bam_stats_samtools", migrate_pytest=True
        )
        subworkflow_create.create()

        with open(subworkflow_dir / "main.nf") as fh:
            new_main_nf = fh.read()
        with open(subworkflow_dir / "meta.yml") as fh:
            new_meta_yml = fh.read()
        nextflow_config = subworkflow_dir / "tests" / "nextflow.config"

        # Check that old files have been copied to the new module
        assert old_main_nf == new_main_nf
        assert old_meta_yml == new_meta_yml
        assert nextflow_config.is_file()

        # Check that pytest folder is deleted
        assert not pytest_dir.is_dir()

        # Check that pytest_modules.yml is updated
        with open(Path(self.nfcore_modules, "tests", "config", "pytest_modules.yml")) as fh:
            modules_yml = yaml.safe_load(fh)
        assert "subworkflows/bam_stats_samtools" not in modules_yml.keys()

    @mock.patch("rich.prompt.Confirm.ask")
    def test_subworkflows_migrate_no_delete(self, mock_rich_ask):
        """Create a subworkflow with the --migrate-pytest option to convert pytest to nf-test.
        Test that pytest directory is not deleted."""
        pytest_dir = Path(self.nfcore_modules, "tests", "subworkflows", "nf-core", "bam_stats_samtools")

        # Clone modules repo with pytests
        shutil.rmtree(self.nfcore_modules)
        Repo.clone_from(GITLAB_URL, self.nfcore_modules, branch=GITLAB_SUBWORKFLOWS_ORG_PATH_BRANCH)

        # Create a module with --migrate-pytest
        mock_rich_ask.return_value = False
        module_create = nf_core.subworkflows.SubworkflowCreate(
            self.nfcore_modules, "bam_stats_samtools", migrate_pytest=True
        )
        module_create.create()

        # Check that pytest folder is not deleted
        assert pytest_dir.is_dir()

        # Check that pytest_modules.yml is updated
        with open(Path(self.nfcore_modules, "tests", "config", "pytest_modules.yml")) as fh:
            modules_yml = yaml.safe_load(fh)
        assert "subworkflows/bam_stats_samtools" not in modules_yml.keys()
