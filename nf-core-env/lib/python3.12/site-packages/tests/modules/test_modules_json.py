import copy
import json
import shutil
from pathlib import Path

from nf_core.components.constants import (
    NF_CORE_MODULES_DEFAULT_BRANCH,
    NF_CORE_MODULES_NAME,
    NF_CORE_MODULES_REMOTE,
)
from nf_core.modules.modules_json import ModulesJson
from nf_core.modules.modules_repo import ModulesRepo
from nf_core.modules.patch import ModulePatch

from ..test_modules import TestModules


class TestModulesCreate(TestModules):
    def test_get_modules_json(self):
        """Checks that the get_modules_json function returns the correct result"""
        mod_json_path = Path(self.pipeline_dir, "modules.json")
        with open(mod_json_path) as fh:
            try:
                mod_json_sb = json.load(fh)
            except json.JSONDecodeError as e:
                raise UserWarning(f"Unable to load JSON file '{mod_json_path}' due to error {e}")

        mod_json_obj = ModulesJson(self.pipeline_dir)
        mod_json = mod_json_obj.get_modules_json()

        # Check that the modules.json hasn't changed
        assert mod_json == mod_json_sb

    def test_mod_json_update(self):
        """Checks whether the update function works properly"""
        mod_json_obj = ModulesJson(self.pipeline_dir)
        # Update the modules.json file
        mod_repo_obj = ModulesRepo()
        mod_json_obj.update("modules", mod_repo_obj, "MODULE_NAME", "GIT_SHA", ["modules"], write_file=False)
        mod_json = mod_json_obj.get_modules_json()
        assert "MODULE_NAME" in mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"]["nf-core"]
        assert "git_sha" in mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"]["nf-core"]["MODULE_NAME"]
        assert "GIT_SHA" == mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"]["nf-core"]["MODULE_NAME"]["git_sha"]
        assert (
            NF_CORE_MODULES_DEFAULT_BRANCH
            == mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"]["nf-core"]["MODULE_NAME"]["branch"]
        )
        assert (
            "modules" in mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"]["nf-core"]["MODULE_NAME"]["installed_by"]
        )

    def test_mod_json_create(self):
        """Test creating a modules.json file from scratch"""
        mod_json_path = Path(self.pipeline_dir, "modules.json")
        # Remove the existing modules.json file
        mod_json_path.unlink()

        # Create the new modules.json file
        # (There are no prompts as long as there are only nf-core modules)
        ModulesJson(self.pipeline_dir).create()

        # Check that the file exists
        assert (mod_json_path).exists()

        # Get the contents of the file
        mod_json_obj = ModulesJson(self.pipeline_dir)
        mod_json = mod_json_obj.get_modules_json()

        mods = ["fastqc", "multiqc"]
        for mod in mods:
            assert mod in mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"]["nf-core"]
            assert "git_sha" in mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"]["nf-core"][mod]
            assert "branch" in mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"]["nf-core"][mod]

    def _modify_main_nf(self, path):
        """Modify a file to test patch creation"""
        with open(path) as fh:
            lines = fh.readlines()
        # Modify $meta.id to $meta.single_end
        lines[1] = '    tag "$meta.single_end"\n'
        with open(path, "w") as fh:
            fh.writelines(lines)

    def test_mod_json_create_with_patch(self):
        """Test creating a modules.json file from scratch when there are patched modules"""
        mod_json_path = Path(self.pipeline_dir, "modules.json")

        # Modify the module
        module_path = Path(self.pipeline_dir, "modules", "nf-core", "fastqc")
        self._modify_main_nf(module_path / "main.nf")

        # Try creating a patch file
        patch_obj = ModulePatch(self.pipeline_dir, NF_CORE_MODULES_REMOTE, NF_CORE_MODULES_DEFAULT_BRANCH)
        patch_obj.patch("fastqc")

        # Remove the existing modules.json file
        mod_json_path.unlink()

        # Create the new modules.json file
        ModulesJson(self.pipeline_dir).create()

        # Check that the file exists
        assert mod_json_path.is_file()

        # Get the contents of the file
        mod_json_obj = ModulesJson(self.pipeline_dir)
        mod_json = mod_json_obj.get_modules_json()

        # Check that fastqc is in the file
        assert "fastqc" in mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"]["nf-core"]
        assert "git_sha" in mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"]["nf-core"]["fastqc"]
        assert "branch" in mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"]["nf-core"]["fastqc"]

        # Check that fastqc/main.nf maintains the changes
        with open(module_path / "main.nf") as fh:
            lines = fh.readlines()
        assert lines[1] == '    tag "$meta.single_end"\n'

    def test_mod_json_up_to_date(self):
        """
        Checks if the modules.json file is up to date
        when no changes have been made to the pipeline
        """
        mod_json_obj = ModulesJson(self.pipeline_dir)
        mod_json_before = mod_json_obj.get_modules_json()
        mod_json_obj.check_up_to_date()
        mod_json_after = mod_json_obj.get_modules_json()

        # Check that the modules.json hasn't changed
        assert mod_json_before == mod_json_after

    def test_mod_json_up_to_date_module_removed(self):
        """
        Reinstall a module that has an entry in the modules.json
        but is missing in the pipeline
        """
        # Remove the fastqc module
        fastqc_path = Path(self.pipeline_dir, "modules", NF_CORE_MODULES_NAME, "fastqc")
        shutil.rmtree(fastqc_path)

        # Check that the modules.json file is up to date, and reinstall the module
        mod_json_obj = ModulesJson(self.pipeline_dir)
        mod_json_obj.check_up_to_date()

        # Check that the module has been reinstalled
        files = ["main.nf", "meta.yml"]
        assert fastqc_path.exists()
        for f in files:
            assert Path(fastqc_path, f).exists()

    def test_mod_json_up_to_date_reinstall_fails(self):
        """
        Try reinstalling a module where the git_sha is invalid
        """
        mod_json_obj = ModulesJson(self.pipeline_dir)

        # Update the fastqc module entry to an invalid git_sha
        mod_json_obj.update("modules", ModulesRepo(), "fastqc", "INVALID_GIT_SHA", ["modules"], write_file=True)

        # Remove the fastqc module
        fastqc_path = Path(self.pipeline_dir, "modules", NF_CORE_MODULES_NAME, "fastqc")
        shutil.rmtree(fastqc_path)

        # Check that the modules.json file is up to date, and remove the fastqc module entry
        mod_json_obj.check_up_to_date()
        mod_json = mod_json_obj.get_modules_json()

        # Check that the module has been removed from the modules.json
        assert "fastqc" not in mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"]["nf-core"]

    def test_mod_json_repo_present(self):
        """Tests the repo_present function"""
        mod_json_obj = ModulesJson(self.pipeline_dir)

        assert mod_json_obj.repo_present(NF_CORE_MODULES_REMOTE) is True
        assert mod_json_obj.repo_present("INVALID_REPO") is False

    def test_mod_json_component_present(self):
        """Tests the component_present function"""
        mod_json_obj = ModulesJson(self.pipeline_dir)

        assert mod_json_obj.component_present("fastqc", NF_CORE_MODULES_REMOTE, NF_CORE_MODULES_NAME, "modules") is True
        assert (
            mod_json_obj.component_present("INVALID_MODULE", NF_CORE_MODULES_REMOTE, NF_CORE_MODULES_NAME, "modules")
            is False
        )
        assert mod_json_obj.component_present("fastqc", "INVALID_REPO", "INVALID_DIR", "modules") is False
        assert mod_json_obj.component_present("INVALID_MODULE", "INVALID_REPO", "INVALID_DIR", "modules") is False

    def test_mod_json_get_module_version(self):
        """Test the get_module_version function"""
        mod_json_obj = ModulesJson(self.pipeline_dir)
        mod_json = mod_json_obj.get_modules_json()
        assert (
            mod_json_obj.get_module_version("fastqc", NF_CORE_MODULES_REMOTE, NF_CORE_MODULES_NAME)
            == mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"]["nf-core"]["fastqc"]["git_sha"]
        )
        assert mod_json_obj.get_module_version("INVALID_MODULE", NF_CORE_MODULES_REMOTE, NF_CORE_MODULES_NAME) is None

    def test_mod_json_dump(self):
        """Tests the dump function"""
        mod_json_obj = ModulesJson(self.pipeline_dir)
        mod_json = mod_json_obj.get_modules_json()
        # Remove the modules.json file
        mod_json_path = Path(self.pipeline_dir, "modules.json")
        mod_json_path.unlink()

        # Check that the dump function creates the file
        mod_json_obj.dump()
        assert mod_json_path.exists()

        # Check that the dump function writes the correct content
        with open(mod_json_path) as f:
            try:
                mod_json_new = json.load(f)
            except json.JSONDecodeError as e:
                raise UserWarning(f"Unable to load JSON file '{mod_json_path}' due to error {e}")
        assert mod_json == mod_json_new

    def test_mod_json_with_empty_modules_value(self):
        # Load module.json and remove the modules entry
        mod_json_obj = ModulesJson(self.pipeline_dir)
        mod_json_obj.create()  # Create modules.json explicitly to get correct module sha
        mod_json_orig = mod_json_obj.get_modules_json()
        mod_json = copy.deepcopy(mod_json_orig)
        mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"] = {}
        # save the altered module.json and load it again to check if it will fix itself
        mod_json_obj.modules_json = mod_json
        mod_json_obj.dump()
        mod_json_obj_new = ModulesJson(self.pipeline_dir)
        mod_json_obj_new.check_up_to_date()
        mod_json_new = mod_json_obj_new.get_modules_json()
        assert mod_json_orig == mod_json_new

    def test_mod_json_with_missing_modules_entry(self):
        # Load module.json and remove the modules entry
        mod_json_obj = ModulesJson(self.pipeline_dir)
        mod_json_obj.create()  # Create modules.json explicitly to get correct module sha
        mod_json_orig = mod_json_obj.get_modules_json()
        mod_json = copy.deepcopy(mod_json_orig)
        mod_json["repos"][NF_CORE_MODULES_REMOTE].pop("modules")
        # save the altered module.json and load it again to check if it will fix itself
        mod_json_obj.modules_json = mod_json
        mod_json_obj.dump()
        mod_json_obj_new = ModulesJson(self.pipeline_dir)
        mod_json_obj_new.check_up_to_date()
        mod_json_new = mod_json_obj_new.get_modules_json()
        assert mod_json_orig == mod_json_new
