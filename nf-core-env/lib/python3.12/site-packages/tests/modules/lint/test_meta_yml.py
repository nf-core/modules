from pathlib import Path

import yaml

import nf_core.modules.lint

from ...test_modules import TestModules


class TestMetaYml(TestModules):
    """Test meta.yml functionality"""

    def test_modules_lint_update_meta_yml(self):
        """update the meta.yml of a module"""
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules, fix=True)
        module_lint.lint(print_results=False, module="bpipe/test")
        assert len(module_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0

    def test_modules_meta_yml_incorrect_licence_field(self):
        """Test linting a module with an incorrect Licence field in meta.yml"""
        with open(self.bpipe_test_module_path / "meta.yml") as fh:
            meta_yml = yaml.safe_load(fh)
        meta_yml["tools"][0]["bpipe"]["licence"] = "[MIT]"
        with open(
            self.bpipe_test_module_path / "meta.yml",
            "w",
        ) as fh:
            fh.write(yaml.dump(meta_yml))
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")

        assert len(module_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) >= 0
        assert len(module_lint.warned) >= 0
        assert module_lint.failed[0].lint_test == "meta_yml_valid"

    def test_modules_meta_yml_output_mismatch(self):
        """Test linting a module with an extra entry in output fields in meta.yml compared to module.output"""
        with open(Path(self.nfcore_modules, "modules", "nf-core", "bpipe", "test", "main.nf")) as fh:
            main_nf = fh.read()
        main_nf_new = main_nf.replace("emit: sequence_report", "emit: bai")
        with open(Path(self.nfcore_modules, "modules", "nf-core", "bpipe", "test", "main.nf"), "w") as fh:
            fh.write(main_nf_new)
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")
        with open(Path(self.nfcore_modules, "modules", "nf-core", "bpipe", "test", "main.nf"), "w") as fh:
            fh.write(main_nf)
        assert len(module_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) >= 0
        assert "Module `meta.yml` does not match `main.nf`" in module_lint.failed[0].message

    def test_modules_meta_yml_incorrect_name(self):
        """Test linting a module with an incorrect name in meta.yml"""
        with open(Path(self.nfcore_modules, "modules", "nf-core", "bpipe", "test", "meta.yml")) as fh:
            meta_yml = yaml.safe_load(fh)
        meta_yml["name"] = "bpipe/test"
        with open(
            Path(self.nfcore_modules, "modules", "nf-core", "bpipe", "test", "meta.yml"),
            "w",
        ) as fh:
            fh.write(yaml.dump(meta_yml))
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")

        assert len(module_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) >= 0
        assert len(module_lint.warned) >= 0
        assert module_lint.failed[0].lint_test == "meta_name"
