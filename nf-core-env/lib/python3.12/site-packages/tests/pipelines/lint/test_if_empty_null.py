from pathlib import Path

import yaml

import nf_core.pipelines.create
import nf_core.pipelines.lint

from ..test_lint import TestLint


class TestLintIfEmptyNull(TestLint):
    def setUp(self) -> None:
        super().setUp()
        self.new_pipeline = self._make_pipeline_copy()
        self.nf_core_yml_path = Path(self.new_pipeline) / ".nf-core.yml"
        with open(self.nf_core_yml_path) as f:
            self.nf_core_yml = yaml.safe_load(f)

    def test_if_empty_null_throws_warn(self):
        """Tests finding ifEmpty(null) in file throws warn in linting"""
        # Create a file and add examples that should fail linting
        txt_file = Path(self.new_pipeline) / "docs" / "test.txt"
        with open(txt_file, "w") as f:
            f.writelines(
                [
                    "ifEmpty(null)\n",
                    "ifEmpty (null)\n",
                    "ifEmpty( null )\n",
                    "ifEmpty ( null )\n",
                    ".ifEmpty(null)\n",
                    ". ifEmpty(null)\n",
                    "|ifEmpty(null)\n",
                    "| ifEmpty(null)\n",
                ]
            )
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj._load()
        result = lint_obj.pipeline_if_empty_null()
        assert len(result["warned"]) == 8
