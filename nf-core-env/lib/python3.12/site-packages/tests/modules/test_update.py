import logging
import shutil
import tempfile
from pathlib import Path
from unittest import mock

import questionary
import yaml

import nf_core.utils
from nf_core.components.constants import NF_CORE_MODULES_NAME, NF_CORE_MODULES_REMOTE
from nf_core.modules.install import ModuleInstall
from nf_core.modules.modules_json import ModulesJson
from nf_core.modules.patch import ModulePatch
from nf_core.modules.update import ModuleUpdate

from ..test_modules import TestModules
from ..utils import (
    GITLAB_BRANCH_TEST_BRANCH,
    GITLAB_BRANCH_TEST_NEW_SHA,
    GITLAB_BRANCH_TEST_OLD_SHA,
    GITLAB_DEFAULT_BRANCH,
    GITLAB_REPO,
    GITLAB_URL,
    OLD_TRIMGALORE_BRANCH,
    OLD_TRIMGALORE_SHA,
    cmp_component,
)


class TestModulesInstall(TestModules):
    def test_install_and_update(self):
        """Installs a module in the pipeline and updates it (no change)"""
        self.mods_install.install("trimgalore")
        update_obj = ModuleUpdate(self.pipeline_dir, show_diff=False)

        # Copy the module files and check that they are unaffected by the update
        tmpdir = Path(tempfile.TemporaryDirectory().name)
        trimgalore_tmpdir = tmpdir / "trimgalore"
        trimgalore_path = Path(self.pipeline_dir, "modules", NF_CORE_MODULES_NAME, "trimgalore")
        shutil.copytree(trimgalore_path, trimgalore_tmpdir)

        assert update_obj.update("trimgalore") is True
        assert cmp_component(trimgalore_tmpdir, trimgalore_path) is True

    def test_install_at_hash_and_update(self):
        """Installs an old version of a module in the pipeline and updates it"""
        assert self.mods_install_old.install("trimgalore")
        update_obj = ModuleUpdate(
            self.pipeline_dir, show_diff=False, update_deps=True, remote_url=GITLAB_URL, branch=OLD_TRIMGALORE_BRANCH
        )

        # Copy the module files and check that they are affected by the update
        tmpdir = Path(tempfile.TemporaryDirectory().name)
        trimgalore_tmpdir = tmpdir / "trimgalore"
        trimgalore_path = Path(self.pipeline_dir, "modules", GITLAB_REPO, "trimgalore")
        shutil.copytree(trimgalore_path, trimgalore_tmpdir)

        assert update_obj.update("trimgalore") is True
        assert cmp_component(trimgalore_tmpdir, trimgalore_path) is False

        # Check that the modules.json is correctly updated
        mod_json_obj = ModulesJson(self.pipeline_dir)
        mod_json = mod_json_obj.get_modules_json()
        # Get the up-to-date git_sha for the module from the ModulesRepo object
        correct_git_sha = update_obj.modules_repo.get_latest_component_version("trimgalore", "modules")
        current_git_sha = mod_json["repos"][GITLAB_URL]["modules"][GITLAB_REPO]["trimgalore"]["git_sha"]
        assert correct_git_sha == current_git_sha

    # Mock questionary answer: do not update module, only show diffs
    @mock.patch.object(questionary.Question, "unsafe_ask", return_value=True)
    def test_install_at_hash_and_update_limit_output(self, mock_prompt):
        """Installs an old version of a module in the pipeline and updates it with limited output reporting"""
        self.caplog.set_level(logging.INFO)
        assert self.mods_install_old.install("trimgalore")

        update_obj = ModuleUpdate(
            self.pipeline_dir,
            show_diff=True,
            update_deps=True,
            remote_url=GITLAB_URL,
            branch=OLD_TRIMGALORE_BRANCH,
            limit_output=True,
        )
        assert update_obj.update("trimgalore")

        # Check changes not shown for non-.nf files
        assert "Changes in 'trimgalore/meta.yml' but not shown" in self.caplog.text
        # Check changes shown for .nf files
        assert "Changes in 'trimgalore/main.nf'" in self.caplog.text
        for line in self.caplog.text.split("\n"):
            if line.startswith("---"):
                assert line.endswith("main.nf")

    def test_install_at_hash_and_update_and_save_diff_to_file(self):
        """Installs an old version of a module in the pipeline and updates it"""
        self.mods_install_old.install("trimgalore")
        patch_path = Path(self.pipeline_dir, "trimgalore.patch")
        update_obj = ModuleUpdate(
            self.pipeline_dir,
            save_diff_fn=patch_path,
            sha=OLD_TRIMGALORE_SHA,
            remote_url=GITLAB_URL,
            branch=OLD_TRIMGALORE_BRANCH,
        )

        # Copy the module files and check that they are affected by the update
        tmpdir = Path(tempfile.TemporaryDirectory().name)
        trimgalore_tmpdir = tmpdir / "trimgalore"
        trimgalore_path = Path(self.pipeline_dir, "modules", GITLAB_REPO, "trimgalore")
        shutil.copytree(trimgalore_path, trimgalore_tmpdir)

        assert update_obj.update("trimgalore") is True
        assert cmp_component(trimgalore_tmpdir, trimgalore_path) is True

        # TODO: Apply the patch to the module

    def test_install_at_hash_and_update_and_save_diff_to_file_limit_output(self):
        """Installs an old version of a module in the pipeline and updates it"""
        # Install old version of trimgalore
        self.mods_install_old.install("trimgalore")
        patch_path = Path(self.pipeline_dir, "trimgalore.patch")
        # Update saving the differences to a patch file and with `limit_output`
        update_obj = ModuleUpdate(
            self.pipeline_dir,
            save_diff_fn=patch_path,
            remote_url=GITLAB_URL,
            branch=OLD_TRIMGALORE_BRANCH,
            limit_output=True,
        )
        assert update_obj.update("trimgalore")

        # Check that the patch file was created
        assert patch_path.exists(), f"Patch file was not created at {patch_path}"

        # Read the contents of the patch file
        with open(patch_path) as fh:
            patch_content = fh.read()
            # Check changes not shown for non-.nf files
            assert "Changes in 'trimgalore/meta.yml' but not shown" in patch_content
            # Check changes only shown for main.nf
            assert "Changes in 'trimgalore/main.nf'" in patch_content
            for line in patch_content:
                if line.startswith("---"):
                    assert line.endswith("main.nf")

    def test_update_all(self):
        """Updates all modules present in the pipeline"""
        update_obj = ModuleUpdate(self.pipeline_dir, update_all=True, show_diff=False)
        # Get the current modules.json
        assert update_obj.update() is True

        # We must reload the modules.json to get the updated version
        mod_json_obj = ModulesJson(self.pipeline_dir)
        mod_json = mod_json_obj.get_modules_json()
        # Loop through all modules and check that they are updated (according to the modules.json file)
        for mod in mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"][NF_CORE_MODULES_NAME]:
            correct_git_sha = list(update_obj.modules_repo.get_component_git_log(mod, "modules", depth=1))[0]["git_sha"]
            current_git_sha = mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"][NF_CORE_MODULES_NAME][mod]["git_sha"]
            assert correct_git_sha == current_git_sha

    def test_update_with_config_fixed_version(self):
        """Try updating when there are entries in the .nf-core.yml"""
        # Install trimgalore at the latest version
        assert self.mods_install_trimgalore.install("trimgalore")

        # Fix the trimgalore version in the .nf-core.yml to an old version
        update_config = {GITLAB_URL: {GITLAB_REPO: {"trimgalore": OLD_TRIMGALORE_SHA}}}
        config_fn, tools_config = nf_core.utils.load_tools_config(self.pipeline_dir)
        setattr(tools_config, "update", update_config)
        assert config_fn is not None and tools_config is not None  # mypy
        with open(Path(self.pipeline_dir, config_fn), "w") as f:
            yaml.dump(tools_config.model_dump(), f)

        # Update all modules in the pipeline
        update_obj = ModuleUpdate(
            self.pipeline_dir, update_all=True, show_diff=False, remote_url=GITLAB_URL, branch=OLD_TRIMGALORE_BRANCH
        )
        assert update_obj.update() is True

        # Check that the git sha for trimgalore is correctly downgraded
        mod_json = ModulesJson(self.pipeline_dir).get_modules_json()
        assert "trimgalore" in mod_json["repos"][GITLAB_URL]["modules"][GITLAB_REPO]
        assert "git_sha" in mod_json["repos"][GITLAB_URL]["modules"][GITLAB_REPO]["trimgalore"]
        assert mod_json["repos"][GITLAB_URL]["modules"][GITLAB_REPO]["trimgalore"]["git_sha"] == OLD_TRIMGALORE_SHA

    def test_update_with_config_dont_update(self):
        """Try updating when module is to be ignored"""
        # Install an old version of trimgalore
        self.mods_install_old.install("trimgalore")

        # Set the trimgalore field to no update in the .nf-core.yml
        update_config = {GITLAB_URL: {GITLAB_REPO: {"trimgalore": False}}}
        config_fn, tools_config = nf_core.utils.load_tools_config(self.pipeline_dir)
        setattr(tools_config, "update", update_config)
        assert config_fn is not None and tools_config is not None  # mypy
        with open(Path(self.pipeline_dir, config_fn), "w") as f:
            yaml.dump(tools_config.model_dump(), f)

        # Update all modules in the pipeline
        update_obj = ModuleUpdate(
            self.pipeline_dir,
            update_all=True,
            show_diff=False,
            sha=OLD_TRIMGALORE_SHA,
            remote_url=GITLAB_URL,
            branch=OLD_TRIMGALORE_BRANCH,
        )
        assert update_obj.update() is True

        # Check that the git sha for trimgalore is correctly downgraded
        mod_json = ModulesJson(self.pipeline_dir).get_modules_json()
        assert "trimgalore" in mod_json["repos"][GITLAB_URL]["modules"][GITLAB_REPO]
        assert "git_sha" in mod_json["repos"][GITLAB_URL]["modules"][GITLAB_REPO]["trimgalore"]
        assert mod_json["repos"][GITLAB_URL]["modules"][GITLAB_REPO]["trimgalore"]["git_sha"] == OLD_TRIMGALORE_SHA

    def test_update_with_config_fix_all(self):
        """Fix the version of all nf-core modules"""
        self.mods_install_trimgalore.install("trimgalore")

        # Fix the version of all nf-core modules in the .nf-core.yml to an old version
        update_config = {GITLAB_URL: OLD_TRIMGALORE_SHA}
        config_fn, tools_config = nf_core.utils.load_tools_config(self.pipeline_dir)
        setattr(tools_config, "update", update_config)
        assert config_fn is not None and tools_config is not None  # mypy
        with open(Path(self.pipeline_dir, config_fn), "w") as f:
            yaml.dump(tools_config.model_dump(), f)

        # Update all modules in the pipeline
        update_obj = ModuleUpdate(
            self.pipeline_dir, update_all=True, show_diff=False, remote_url=GITLAB_URL, branch=OLD_TRIMGALORE_BRANCH
        )
        assert update_obj.update() is True

        # Check that the git sha for trimgalore is correctly downgraded
        mod_json = ModulesJson(self.pipeline_dir).get_modules_json()
        assert "git_sha" in mod_json["repos"][GITLAB_URL]["modules"][GITLAB_REPO]["trimgalore"]
        assert mod_json["repos"][GITLAB_URL]["modules"][GITLAB_REPO]["trimgalore"]["git_sha"] == OLD_TRIMGALORE_SHA

    def test_update_with_config_no_updates(self):
        """Don't update any nf-core modules"""
        assert self.mods_install_old.install("trimgalore")
        old_mod_json = ModulesJson(self.pipeline_dir).get_modules_json()

        # Fix the version of all nf-core modules in the .nf-core.yml to an old version
        update_config = {GITLAB_URL: False}
        config_fn, tools_config = nf_core.utils.load_tools_config(self.pipeline_dir)
        setattr(tools_config, "update", update_config)
        assert config_fn is not None and tools_config is not None  # mypy
        with open(Path(self.pipeline_dir, config_fn), "w") as f:
            yaml.dump(tools_config.model_dump(), f)

        # Update all modules in the pipeline
        update_obj = ModuleUpdate(
            self.pipeline_dir,
            update_all=True,
            show_diff=False,
            sha=OLD_TRIMGALORE_SHA,
            remote_url=GITLAB_URL,
            branch=OLD_TRIMGALORE_BRANCH,
        )
        assert update_obj.update() is True

        # Check that the git sha for trimgalore is correctly downgraded and none of the modules has changed
        mod_json = ModulesJson(self.pipeline_dir).get_modules_json()
        for module in mod_json["repos"][GITLAB_URL]["modules"][GITLAB_REPO]:
            assert "git_sha" in mod_json["repos"][GITLAB_URL]["modules"][GITLAB_REPO][module]
            assert (
                mod_json["repos"][GITLAB_URL]["modules"][GITLAB_REPO][module]["git_sha"]
                == old_mod_json["repos"][GITLAB_URL]["modules"][GITLAB_REPO][module]["git_sha"]
            )

    def test_update_different_branch_single_module(self):
        """Try updating a module in a specific branch"""
        install_obj = ModuleInstall(
            self.pipeline_dir,
            prompt=False,
            force=False,
            remote_url=GITLAB_URL,
            branch=GITLAB_BRANCH_TEST_BRANCH,
            sha=GITLAB_BRANCH_TEST_OLD_SHA,
        )
        assert install_obj.install("fastp")

        update_obj = ModuleUpdate(
            self.pipeline_dir,
            update_deps=True,
            remote_url=GITLAB_URL,
            branch=GITLAB_BRANCH_TEST_BRANCH,
            show_diff=False,
        )
        update_obj.update("fastp")

        # Verify that the branch entry was updated correctly
        modules_json = ModulesJson(self.pipeline_dir)
        assert (
            modules_json.get_component_branch(self.component_type, "fastp", GITLAB_URL, GITLAB_REPO)
            == GITLAB_BRANCH_TEST_BRANCH
        )
        assert modules_json.get_module_version("fastp", GITLAB_URL, GITLAB_REPO) == GITLAB_BRANCH_TEST_NEW_SHA

    def test_update_different_branch_mixed_modules_main(self):
        """Try updating all modules where MultiQC is installed from main branch"""
        # Install fastp
        assert self.mods_install_gitlab_old.install("fastp")

        # Install MultiQC from gitlab default branch
        assert self.mods_install_gitlab.install("multiqc")

        # Try updating
        update_obj = ModuleUpdate(self.pipeline_dir, update_all=True, show_diff=False)
        assert update_obj.update() is True

        modules_json = ModulesJson(self.pipeline_dir)
        # Verify that the branch entry was updated correctly
        assert (
            modules_json.get_component_branch(self.component_type, "fastp", GITLAB_URL, GITLAB_REPO)
            == GITLAB_BRANCH_TEST_BRANCH
        )
        assert modules_json.get_module_version("fastp", GITLAB_URL, GITLAB_REPO) == GITLAB_BRANCH_TEST_NEW_SHA
        # MultiQC is present in both branches but should've been updated using the 'main' branch
        assert (
            modules_json.get_component_branch(self.component_type, "multiqc", GITLAB_URL, GITLAB_REPO)
            == GITLAB_DEFAULT_BRANCH
        )

    def test_update_different_branch_mix_modules_branch_test(self):
        """Try updating all modules where MultiQC is installed from branch-test branch"""
        # Install multiqc from the branch-test branch
        assert self.mods_install_gitlab_old.install(
            "multiqc"
        )  # Force as the same module is installed from github nf-core modules repo
        modules_json = ModulesJson(self.pipeline_dir)
        update_obj = ModuleUpdate(
            self.pipeline_dir,
            update_all=True,
            show_diff=False,
            remote_url=GITLAB_URL,
            branch=GITLAB_BRANCH_TEST_BRANCH,
            sha=GITLAB_BRANCH_TEST_NEW_SHA,
        )
        assert update_obj.update()

        assert (
            modules_json.get_component_branch(self.component_type, "multiqc", GITLAB_URL, GITLAB_REPO)
            == GITLAB_BRANCH_TEST_BRANCH
        )
        assert modules_json.get_module_version("multiqc", GITLAB_URL, GITLAB_REPO) == GITLAB_BRANCH_TEST_NEW_SHA

    # Mock questionary answer: do not update module, only show diffs
    @mock.patch.object(questionary.Question, "unsafe_ask", return_value=False)
    def test_update_only_show_differences(self, mock_prompt):
        """Try updating all modules showing differences.
        Only show diffs, don't actually save any updated files.
        Check that the sha in modules.json is not changed."""

        # Update modules to a fixed old SHA
        update_old = ModuleUpdate(
            self.pipeline_dir, update_all=True, show_diff=False, sha="5e34754d42cd2d5d248ca8673c0a53cdf5624905"
        )
        update_old.update()

        tmpdir = Path(tempfile.TemporaryDirectory().name)
        shutil.copytree(Path(self.pipeline_dir, "modules", NF_CORE_MODULES_NAME), tmpdir)

        update_obj = ModuleUpdate(self.pipeline_dir, update_all=True, show_diff=True)
        assert ModuleUpdate(self.pipeline_dir, update_all=True, show_diff=True).update()

        mod_json = ModulesJson(self.pipeline_dir).get_modules_json()
        # Loop through all modules and check that they are NOT updated (according to the modules.json file)
        # A module that can be updated but shouldn't is fastqc
        # Module multiqc is already up to date so don't check
        mod = "fastqc"
        non_updated_git_sha = list(update_obj.modules_repo.get_component_git_log(mod, "modules", depth=1))[0]["git_sha"]
        current_git_sha = mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"][NF_CORE_MODULES_NAME][mod]["git_sha"]
        assert non_updated_git_sha != current_git_sha
        assert cmp_component(Path(tmpdir, mod), Path(self.pipeline_dir, "modules", NF_CORE_MODULES_NAME, mod)) is True

    # Mock questionary answer: do not update module, only show diffs
    @mock.patch.object(questionary.Question, "unsafe_ask", return_value=False)
    def test_update_only_show_differences_when_patch(self, mock_prompt):
        """Try updating all modules showing differences when there's a patched module.
        Don't update some of them.
        Check that the sha in modules.json is not changed."""
        modules_json = ModulesJson(self.pipeline_dir)
        update_obj = ModuleUpdate(self.pipeline_dir, update_all=True, show_diff=True)

        # Update modules to a fixed old SHA
        update_old = ModuleUpdate(
            self.pipeline_dir, update_all=True, show_diff=False, sha="5e34754d42cd2d5d248ca8673c0a53cdf5624905"
        )
        assert update_old.update()

        # Modify fastqc module, it will have a patch which will be applied during update
        # We modify fastqc because it's one of the modules that can be updated and there's another one before it (custom/dumpsoftwareversions)
        module_path = Path(self.pipeline_dir, "modules", "nf-core", "fastqc")
        main_path = Path(module_path, "main.nf")
        with open(main_path) as fh:
            lines = fh.readlines()
        for line_index in range(len(lines)):
            if lines[line_index] == "    label 'process_medium'\n":
                lines[line_index] = "    label 'process_low'\n"
                break
        with open(main_path, "w") as fh:
            fh.writelines(lines)
        # Create a patch file
        patch_obj = ModulePatch(self.pipeline_dir)
        patch_obj.patch("fastqc")
        # Check that a patch file with the correct name has been created
        assert "fastqc.diff" in [f.name for f in module_path.glob("*.diff")]

        # Update all modules
        assert update_obj.update() is True

        mod_json = modules_json.get_modules_json()
        # Loop through all modules and check that they are NOT updated (according to the modules.json file)
        # A module that can be updated but shouldn't is fastqc
        # Module multiqc is already up to date so don't check
        mod = "fastqc"
        correct_git_sha = list(update_obj.modules_repo.get_component_git_log(mod, "modules", depth=1))[0]["git_sha"]
        current_git_sha = mod_json["repos"][NF_CORE_MODULES_REMOTE]["modules"][NF_CORE_MODULES_NAME][mod]["git_sha"]
        assert correct_git_sha != current_git_sha

    def test_update_module_with_extra_config_file(self):
        """Try updating a module with a config file"""
        # Install the module
        assert self.mods_install.install("trimgalore")
        # Add a nextflow_test.config file to the module
        trimgalore_path = Path(self.pipeline_dir, "modules", "nf-core", "trimgalore")
        Path(trimgalore_path, "nextflow_test.config").touch()
        with open(Path(trimgalore_path, "nextflow_test.config"), "w") as fh:
            fh.write("params.my_param = 'my_value'\n")
        # Update the module
        update_obj = ModuleUpdate(self.pipeline_dir, show_diff=False)
        assert update_obj.update("trimgalore")
        # Check that the nextflow_test.config file is still there
        assert Path(trimgalore_path, "nextflow_test.config").exists()
        with open(Path(trimgalore_path, "nextflow_test.config")) as fh:
            assert "params.my_param = 'my_value'" in fh.read()
