from pathlib import Path

import nf_core.pipelines.lint

from ..test_lint import TestLint


class TestLintLocalComponentStructure(TestLint):
    def setUp(self) -> None:
        super().setUp()
        self.new_pipeline = self._make_pipeline_copy()

    def test_local_component_structure(self):
        local_modules = Path(self.new_pipeline, "modules", "local")
        local_swf = Path(self.new_pipeline, "subworkflows", "local")
        local_modules.mkdir(parents=True, exist_ok=True)
        local_swf.mkdir(parents=True, exist_ok=True)

        (local_modules / "dummy_module.nf").touch()
        (local_swf / "dummy_subworkflow.nf").touch()

        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj._load()

        results = lint_obj.local_component_structure()
        assert len(results.get("warned", [])) == 2
        assert len(results.get("failed", [])) == 0
        assert len(results.get("ignored", [])) == 0
