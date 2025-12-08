"""Tests covering the modules commands"""

import json
import unittest
from pathlib import Path

import pytest
import requests_cache
import responses
import ruamel.yaml

import nf_core.modules
import nf_core.modules.create
import nf_core.modules.install
import nf_core.modules.modules_repo
import nf_core.modules.remove
import nf_core.pipelines.create.create
from nf_core import __version__
from nf_core.pipelines.lint_utils import run_prettier_on_file
from nf_core.utils import NFCoreYamlConfig

from .utils import (
    GITLAB_BRANCH_TEST_BRANCH,
    GITLAB_BRANCH_TEST_OLD_SHA,
    GITLAB_DEFAULT_BRANCH,
    GITLAB_URL,
    OLD_TRIMGALORE_BRANCH,
    OLD_TRIMGALORE_SHA,
    create_tmp_pipeline,
    mock_anaconda_api_calls,
    mock_biocontainers_api_calls,
    mock_biotools_api_calls,
)


def create_modules_repo_dummy(tmp_dir):
    """Create a dummy copy of the nf-core/modules repo"""
    yaml = ruamel.yaml.YAML()
    yaml.preserve_quotes = True
    yaml.indent(mapping=2, sequence=2, offset=0)

    root_dir = Path(tmp_dir, "modules")
    Path(root_dir, "modules", "nf-core").mkdir(parents=True)
    Path(root_dir, "tests", "modules", "nf-core").mkdir(parents=True)
    Path(root_dir, "tests", "config").mkdir(parents=True)

    nf_core_yml = NFCoreYamlConfig(nf_core_version=__version__, repository_type="modules", org_path="nf-core")
    with open(Path(root_dir, ".nf-core.yml"), "w") as fh:
        yaml.dump(nf_core_yml.model_dump(), fh)
    # mock biocontainers and anaconda response and biotools response
    with responses.RequestsMock() as rsps:
        mock_anaconda_api_calls(rsps, "bpipe", "0.9.13--hdfd78af_0")
        mock_biocontainers_api_calls(rsps, "bpipe", "0.9.13--hdfd78af_0")
        mock_biotools_api_calls(rsps, "bpipe")
        # bpipe is a valid package on bioconda that is very unlikely to ever be added to nf-core/modules
        module_create = nf_core.modules.create.ModuleCreate(
            root_dir, "bpipe/test", "@author", "process_single", True, False
        )
        with requests_cache.disabled():
            assert module_create.create()

    # Remove doi from meta.yml which makes lint fail
    meta_yml_path = Path(root_dir, "modules", "nf-core", "bpipe", "test", "meta.yml")

    with open(str(meta_yml_path)) as fh:
        meta_yml = yaml.load(fh)
    del meta_yml["tools"][0]["bpipe"]["doi"]
    meta_yml["keywords"] = ["pipelines", "bioinformatics", "run"]
    with open(str(meta_yml_path), "w") as fh:
        yaml.dump(meta_yml, fh)
        run_prettier_on_file(fh.name)
    # Add dummy content to main.nf.test.snap
    test_snap_path = Path(root_dir, "modules", "nf-core", "bpipe", "test", "tests", "main.nf.test.snap")

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

    # remove "TODO" statements from main.nf
    main_nf_path = Path(root_dir, "modules", "nf-core", "bpipe", "test", "main.nf")
    with open(main_nf_path) as fh:
        main_nf = fh.read()
    main_nf = main_nf.replace("TODO", "")
    with open(main_nf_path, "w") as fh:
        fh.write(main_nf)

    # remove "TODO" statements from main.nf.test
    main_nf_test_path = Path(root_dir, "modules", "nf-core", "bpipe", "test", "tests", "main.nf.test")
    with open(main_nf_test_path) as fh:
        main_nf_test = fh.read()
    main_nf_test = main_nf_test.replace("TODO", "")
    with open(main_nf_test_path, "w") as fh:
        fh.write(main_nf_test)

    return root_dir


class TestModules(unittest.TestCase):
    """Class for modules tests"""

    def setUp(self):
        """Create a new PipelineSchema and Launch objects"""
        self.component_type = "modules"

        # Set up the schema
        self.tmp_dir, self.template_dir, self.pipeline_name, self.pipeline_dir = create_tmp_pipeline()

        # Set up install objects
        self.mods_install = nf_core.modules.install.ModuleInstall(self.pipeline_dir, prompt=False, force=True)
        self.mods_install_old = nf_core.modules.install.ModuleInstall(
            self.pipeline_dir,
            prompt=False,
            force=False,
            sha=OLD_TRIMGALORE_SHA,
            remote_url=GITLAB_URL,
            branch=OLD_TRIMGALORE_BRANCH,
        )
        self.mods_install_trimgalore = nf_core.modules.install.ModuleInstall(
            self.pipeline_dir,
            prompt=False,
            force=False,
            remote_url=GITLAB_URL,
            branch=OLD_TRIMGALORE_BRANCH,
        )
        self.mods_install_gitlab = nf_core.modules.install.ModuleInstall(
            self.pipeline_dir,
            prompt=False,
            force=False,
            remote_url=GITLAB_URL,
            branch=GITLAB_DEFAULT_BRANCH,
        )
        self.mods_install_gitlab_old = nf_core.modules.install.ModuleInstall(
            self.pipeline_dir,
            prompt=False,
            force=False,
            remote_url=GITLAB_URL,
            branch=GITLAB_BRANCH_TEST_BRANCH,
            sha=GITLAB_BRANCH_TEST_OLD_SHA,
        )

        # Set up remove objects
        self.mods_remove = nf_core.modules.remove.ModuleRemove(self.pipeline_dir)
        self.mods_remove_gitlab = nf_core.modules.remove.ModuleRemove(
            self.pipeline_dir,
            remote_url=GITLAB_URL,
            branch=GITLAB_DEFAULT_BRANCH,
        )

        # Set up the nf-core/modules repo dummy
        self.nfcore_modules = create_modules_repo_dummy(self.tmp_dir)

        # Common path to the bpipe/test module used in tests
        self.bpipe_test_module_path = Path(self.nfcore_modules, "modules", "nf-core", "bpipe", "test")

    def test_modulesrepo_class(self):
        """Initialise a modules repo object"""
        modrepo = nf_core.modules.modules_repo.ModulesRepo()
        assert modrepo.repo_path == "nf-core"
        assert modrepo.branch == "master"

    @pytest.fixture(autouse=True)
    def _use_caplog(self, caplog):
        self.caplog = caplog
