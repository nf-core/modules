"""Tests covering the subworkflows commands"""

import json
import unittest
from pathlib import Path

import pytest

import nf_core.modules
import nf_core.modules.install
import nf_core.pipelines.create.create
import nf_core.subworkflows

from .utils import (
    CROSS_ORGANIZATION_BRANCH,
    CROSS_ORGANIZATION_URL,
    GITLAB_SUBWORKFLOWS_BRANCH,
    GITLAB_SUBWORKFLOWS_ORG_PATH_BRANCH,
    GITLAB_URL,
    OLD_SUBWORKFLOWS_SHA,
    create_tmp_pipeline,
)


def create_modules_repo_dummy(tmp_dir):
    """Create a dummy copy of the nf-core/modules repo"""

    root_dir = Path(tmp_dir, "modules")
    Path(root_dir, "modules").mkdir(parents=True, exist_ok=True)
    Path(root_dir, "subworkflows", "nf-core").mkdir(parents=True, exist_ok=True)
    Path(root_dir, "tests", "config").mkdir(parents=True, exist_ok=True)
    with open(Path(root_dir, ".nf-core.yml"), "w") as fh:
        fh.writelines(["repository_type: modules", "\n", "org_path: nf-core", "\n"])
    # TODO Add a mock here
    subworkflow_create = nf_core.subworkflows.SubworkflowCreate(root_dir, "test_subworkflow", "@author", True)
    subworkflow_create.create()

    # Add dummy content to main.nf.test.snap
    test_snap_path = Path(
        root_dir,
        "subworkflows",
        "nf-core",
        "test_subworkflow",
        "tests",
        "main.nf.test.snap",
    )
    test_snap_path.parent.mkdir(parents=True, exist_ok=True)
    with open(test_snap_path, "w") as fh:
        json.dump(
            {
                "my test": {
                    "content": [
                        {
                            "0": [],
                            "versions": {},
                        }
                    ]
                }
            },
            fh,
            indent=4,
        )

    return root_dir


class TestSubworkflows(unittest.TestCase):
    """Class for subworkflows tests"""

    def setUp(self):
        """Create a new PipelineStructure and Launch objects"""
        self.component_type = "subworkflows"

        # Set up the pipeline structure
        self.tmp_dir, self.template_dir, self.pipeline_name, self.pipeline_dir = create_tmp_pipeline()

        # Set up the nf-core/modules repo dummy
        self.nfcore_modules = create_modules_repo_dummy(self.tmp_dir)

        # Set up install objects
        self.subworkflow_install = nf_core.subworkflows.SubworkflowInstall(self.pipeline_dir, prompt=False, force=False)
        self.subworkflow_install_gitlab = nf_core.subworkflows.SubworkflowInstall(
            self.pipeline_dir,
            prompt=False,
            force=False,
            remote_url=GITLAB_URL,
            branch=GITLAB_SUBWORKFLOWS_BRANCH,
        )
        self.subworkflow_install_gitlab_same_org_path = nf_core.subworkflows.SubworkflowInstall(
            self.pipeline_dir,
            prompt=False,
            force=False,
            remote_url=GITLAB_URL,
            branch=GITLAB_SUBWORKFLOWS_ORG_PATH_BRANCH,
        )
        self.subworkflow_install_old = nf_core.subworkflows.SubworkflowInstall(
            self.pipeline_dir,
            prompt=False,
            force=False,
            sha=OLD_SUBWORKFLOWS_SHA,
        )
        self.subworkflow_install_module_change = nf_core.subworkflows.SubworkflowInstall(
            self.pipeline_dir,
            prompt=False,
            force=False,
            sha="8c343b3c8a0925949783dc547666007c245c235b",
        )
        self.subworkflow_install_cross_org = nf_core.subworkflows.SubworkflowInstall(
            self.pipeline_dir, remote_url=CROSS_ORGANIZATION_URL, branch=CROSS_ORGANIZATION_BRANCH
        )
        # Another instance to avoid cross-contamination
        self.subworkflow_install_cross_org_again = nf_core.subworkflows.SubworkflowInstall(
            self.pipeline_dir, remote_url=CROSS_ORGANIZATION_URL, branch=CROSS_ORGANIZATION_BRANCH
        )

        self.mods_install = nf_core.modules.install.ModuleInstall(self.pipeline_dir, prompt=False, force=True)

        # Set up remove objects
        self.subworkflow_remove = nf_core.subworkflows.SubworkflowRemove(self.pipeline_dir)
        self.subworkflow_remove_cross_org = nf_core.subworkflows.SubworkflowRemove(
            self.pipeline_dir, remote_url=CROSS_ORGANIZATION_URL, branch=CROSS_ORGANIZATION_BRANCH
        )

    @pytest.fixture(autouse=True)
    def _use_caplog(self, caplog):
        self.caplog = caplog
