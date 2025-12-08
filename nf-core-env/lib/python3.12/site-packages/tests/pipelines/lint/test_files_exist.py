from pathlib import Path

from ruamel.yaml import YAML

import nf_core.pipelines.lint

from ..test_lint import TestLint


class TestLintFilesExist(TestLint):
    def setUp(self) -> None:
        super().setUp()
        self.new_pipeline = self._make_pipeline_copy()
        self.lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)

    def test_files_exist_missing_config(self):
        """Lint test: critical files missing FAIL"""

        Path(self.new_pipeline, "CHANGELOG.md").unlink()

        assert self.lint_obj._load()
        self.lint_obj.nf_config["manifest.name"] = "nf-core/testpipeline"

        results = self.lint_obj.files_exist()
        assert "File not found: `CHANGELOG.md`" in results["failed"]

    def test_files_exist_missing_main(self):
        """Check if missing main issues warning"""

        Path(self.new_pipeline, "main.nf").unlink()

        assert self.lint_obj._load()

        results = self.lint_obj.files_exist()
        assert "File not found: `main.nf`" in results["warned"]

    def test_files_exist_deprecated_file(self):
        """Check whether deprecated file issues warning"""

        Path(self.new_pipeline, "parameters.settings.json").touch()

        assert self.lint_obj._load()

        results = self.lint_obj.files_exist()
        assert results["failed"] == ["File must be removed: `parameters.settings.json`"]

    def test_files_exist_pass(self):
        """Lint check should pass if all files are there"""

        assert self.lint_obj._load()

        results = self.lint_obj.files_exist()
        assert results["failed"] == []

    def test_files_exist_pass_conditional_nfschema(self):
        # replace nf-validation with nf-schema in nextflow.config
        with open(Path(self.new_pipeline, "nextflow.config")) as f:
            config = f.read()
        config = config.replace("nf-validation", "nf-schema")
        with open(Path(self.new_pipeline, "nextflow.config"), "w") as f:
            f.write(config)

        assert self.lint_obj._load()
        self.lint_obj.nf_config["manifest.schema"] = "nf-core"
        results = self.lint_obj.files_exist()
        assert results["failed"] == []
        assert results["ignored"] == []

    def test_files_exists_pass_nf_core_yml_config(self):
        """Check if linting passes with a valid nf-core.yml config"""
        valid_yaml = """
        files_exist:
            - .github/CONTRIBUTING.md
            - CITATIONS.md
        """
        yaml = YAML()
        nf_core_yml_path = Path(self.new_pipeline, ".nf-core.yml")
        nf_core_yml = yaml.load(nf_core_yml_path)

        nf_core_yml["lint"] = yaml.load(valid_yaml)
        yaml.dump(nf_core_yml, nf_core_yml_path)

        self.lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        assert self.lint_obj._load()

        results = self.lint_obj.files_exist()
        assert results["failed"] == []
        assert "File is ignored: `.github/CONTRIBUTING.md`" in results["ignored"]
        assert "File is ignored: `CITATIONS.md`" in results["ignored"]

    def test_files_exists_fail_nf_core_yml_config(self):
        """Check if linting fails with a valid nf-core.yml config"""
        valid_yaml = """
        files_exist:
            - CITATIONS.md
        """

        # remove CITATIONS.md
        Path(self.new_pipeline, "CITATIONS.md").unlink()
        assert self.lint_obj._load()
        # test first if linting fails correctly
        results = self.lint_obj.files_exist()
        assert "File not found: `CITATIONS.md`" in results["failed"]

        yaml = YAML()
        nf_core_yml_path = Path(self.new_pipeline, ".nf-core.yml")
        nf_core_yml = yaml.load(nf_core_yml_path)

        nf_core_yml["lint"] = yaml.load(valid_yaml)
        yaml.dump(nf_core_yml, nf_core_yml_path)

        self.lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        assert self.lint_obj._load()

        results = self.lint_obj.files_exist()
        assert results["failed"] == []
        assert "File is ignored: `CITATIONS.md`" in results["ignored"]
