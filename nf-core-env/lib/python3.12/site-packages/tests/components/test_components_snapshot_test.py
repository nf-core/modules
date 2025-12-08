"""Test the 'modules test' or 'subworkflows test' command which runs nf-test test."""

import shutil
from pathlib import Path

import pytest

from nf_core.components.components_test import ComponentsTest
from nf_core.utils import set_wd

from ..test_components import TestComponents


class TestTestComponentsUtils(TestComponents):
    def test_components_test_check_inputs(self):
        """Test the check_inputs() function - raise UserWarning because module doesn't exist"""
        with set_wd(self.nfcore_modules):
            meta_builder = ComponentsTest(component_type="modules", component_name="none", no_prompts=True)
            with pytest.raises(UserWarning) as excinfo:
                meta_builder.check_inputs()
        assert "Cannot find directory" in str(excinfo.value)

    def test_components_test_no_name_no_prompts(self):
        """Test the check_inputs() function - raise UserWarning prompts are deactivated and module name is not provided."""
        with set_wd(self.nfcore_modules):
            meta_builder = ComponentsTest(component_type="modules", component_name=None, no_prompts=True)
            with pytest.raises(UserWarning) as excinfo:
                meta_builder.check_inputs()
        assert "Module name not provided and prompts deactivated." in str(excinfo.value)

    def test_components_test_no_installed_modules(self):
        """Test the check_inputs() function - raise UserWarning because installed modules were not found"""
        with set_wd(self.nfcore_modules):
            module_dir = Path(self.nfcore_modules, "modules")
            shutil.rmtree(module_dir)
            module_dir.mkdir()
            meta_builder = ComponentsTest(component_type="modules", component_name=None, no_prompts=False)
            meta_builder.repo_type = "modules"
            with pytest.raises(LookupError) as excinfo:
                meta_builder.check_inputs()
        assert "Nothing installed from" in str(excinfo.value)
