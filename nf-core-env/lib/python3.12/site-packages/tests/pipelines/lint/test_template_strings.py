import subprocess
from pathlib import Path

import yaml

import nf_core.pipelines.create
import nf_core.pipelines.lint

from ..test_lint import TestLint


class TestLintTemplateStrings(TestLint):
    def setUp(self) -> None:
        super().setUp()
        self.new_pipeline = self._make_pipeline_copy()
        self.nf_core_yml_path = Path(self.new_pipeline) / ".nf-core.yml"
        with open(self.nf_core_yml_path) as f:
            self.nf_core_yml = yaml.safe_load(f)

    def test_template_strings(self):
        """Tests finding a template string in a file fails linting."""
        # Add template string to a file
        txt_file = Path(self.new_pipeline) / "docs" / "test.txt"
        with open(txt_file, "w") as f:
            f.write("my {{ template_string }}")
        subprocess.check_output(["git", "add", "docs"], cwd=self.new_pipeline)
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj._load()
        result = lint_obj.template_strings()
        assert len(result["failed"]) == 1
        assert len(result["ignored"]) == 0

    def test_template_strings_ignored(self):
        """Tests ignoring template_strings"""
        # Ignore template_strings test
        valid_yaml = """
        template_strings: false
        """
        self.nf_core_yml["lint"] = yaml.safe_load(valid_yaml)
        with open(self.nf_core_yml_path, "w") as f:
            yaml.safe_dump(self.nf_core_yml, f)
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj._load()
        lint_obj._lint_pipeline()
        assert len(lint_obj.failed) == 0
        assert len(lint_obj.ignored) == 1

    def test_template_strings_ignore_file(self):
        """Tests ignoring template_strings file"""
        # Add template string to a file
        txt_file = Path(self.new_pipeline) / "docs" / "test.txt"
        with open(txt_file, "w") as f:
            f.write("my {{ template_string }}")

        subprocess.check_output(["git", "add", "docs"], cwd=self.new_pipeline)

        # Ignore template_strings test
        valid_yaml = """
        template_strings:
            - docs/test.txt
        """
        self.nf_core_yml["lint"] = yaml.safe_load(valid_yaml)
        with open(self.nf_core_yml_path, "w") as f:
            yaml.safe_dump(self.nf_core_yml, f)

        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj._load()
        result = lint_obj.template_strings()

        assert len(result["failed"]) == 0
        assert len(result["ignored"]) == 1
