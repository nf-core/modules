import os
import re
import tempfile
from pathlib import Path

import pytest
import ruamel.yaml

import nf_core.modules.bump_versions
from nf_core import __version__
from nf_core.modules.modules_utils import ModuleExceptionError
from nf_core.utils import NFCoreYamlConfig

from ..test_modules import TestModules


class TestModulesBumpVersions(TestModules):
    def test_modules_bump_versions_single_module(self):
        """Test updating a single module"""
        # Change the bpipe/test version to an older version
        env_yml_path = os.path.join(self.nfcore_modules, "modules", "nf-core", "bpipe", "test", "environment.yml")
        with open(env_yml_path) as fh:
            content = fh.read()
        new_content = re.sub(r"bioconda::star=\d.\d.\d\D?", r"bioconda::star=2.6.1d", content)
        with open(env_yml_path, "w") as fh:
            fh.write(new_content)
        version_bumper = nf_core.modules.bump_versions.ModuleVersionBumper(pipeline_dir=self.nfcore_modules)
        modules = version_bumper.bump_versions(module="bpipe/test")
        assert len(version_bumper.failed) == 0
        assert [m.component_name for m in modules] == ["bpipe/test"]

    def test_modules_bump_versions_all_modules(self):
        """Test updating all modules"""
        version_bumper = nf_core.modules.bump_versions.ModuleVersionBumper(pipeline_dir=self.nfcore_modules)
        modules = version_bumper.bump_versions(all_modules=True)
        assert len(version_bumper.failed) == 0
        assert [m.component_name for m in modules] == ["bpipe/test"]

    @staticmethod
    def _mock_nf_core_yml(root_dir: Path) -> None:
        """Mock the .nf_core.yml"""
        yaml = ruamel.yaml.YAML()
        yaml.preserve_quotes = True
        yaml.indent(mapping=2, sequence=2, offset=0)
        nf_core_yml = NFCoreYamlConfig(nf_core_version=__version__, repository_type="modules", org_path="nf-core")
        with open(Path(root_dir, ".nf-core.yml"), "w") as fh:
            yaml.dump(nf_core_yml.model_dump(), fh)

    @staticmethod
    def _mock_modules(root_dir: Path, modules: list[str]) -> None:
        """Mock the directory for a given module (or sub-module) for use with `dry_run`"""
        nf_core_dir = root_dir / "modules" / "nf-core"
        for module in modules:
            if "/" in module:
                module, sub_module = module.split("/")
                module_dir = nf_core_dir / module / sub_module
            else:
                module_dir = nf_core_dir / module
            module_dir.mkdir(parents=True)
            module_main = module_dir / "main.nf"
            with module_main.open("w"):
                pass

    def test_modules_bump_versions_multiple_modules(self):
        """Test updating all modules when multiple modules are present"""
        # mock the fgbio directory
        root_dir = Path(tempfile.TemporaryDirectory().name)
        self._mock_modules(root_dir=root_dir, modules=["fqgrep", "fqtk"])
        # mock the ".nf-core.yml"
        self._mock_nf_core_yml(root_dir=root_dir)

        # run it with dryrun to return the modules that it found
        version_bumper = nf_core.modules.bump_versions.ModuleVersionBumper(pipeline_dir=root_dir)
        modules = version_bumper.bump_versions(all_modules=True, dry_run=True)
        assert sorted([m.component_name for m in modules]) == sorted(["fqgrep", "fqtk"])

    def test_modules_bump_versions_submodules(self):
        """Test updating a submodules"""
        # mock the fgbio directory
        root_dir = Path(tempfile.TemporaryDirectory().name)
        in_modules = ["fgbio/callduplexconsensusreads", "fgbio/groupreadsbyumi"]
        self._mock_modules(root_dir=root_dir, modules=in_modules)
        # mock the ".nf-core.yml"
        self._mock_nf_core_yml(root_dir=root_dir)

        # run it with dryrun to return the modules that it found
        version_bumper = nf_core.modules.bump_versions.ModuleVersionBumper(pipeline_dir=root_dir)
        out_modules = version_bumper.bump_versions(module="fgbio", dry_run=True)
        assert sorted([m.component_name for m in out_modules]) == sorted(in_modules)

    def test_modules_bump_versions_fail(self):
        """Fail updating a module with wrong name"""
        version_bumper = nf_core.modules.bump_versions.ModuleVersionBumper(pipeline_dir=self.nfcore_modules)
        with pytest.raises(ModuleExceptionError) as excinfo:
            version_bumper.bump_versions(module="no/module")
        assert "Could not find the specified module:" in str(excinfo.value)

    def test_modules_bump_versions_fail_unknown_version(self):
        """Fail because of an unknown version"""
        # Change the bpipe/test version to an older version
        env_yml_path = os.path.join(self.nfcore_modules, "modules", "nf-core", "bpipe", "test", "environment.yml")
        with open(env_yml_path) as fh:
            content = fh.read()
        new_content = re.sub(r"bioconda::bpipe=\d.\d.\d\D?", r"bioconda::bpipe=xxx", content)
        with open(env_yml_path, "w") as fh:
            fh.write(new_content)
        version_bumper = nf_core.modules.bump_versions.ModuleVersionBumper(pipeline_dir=self.nfcore_modules)
        version_bumper.bump_versions(module="bpipe/test")
        assert "Conda package had unknown version" in version_bumper.failed[0][0]
