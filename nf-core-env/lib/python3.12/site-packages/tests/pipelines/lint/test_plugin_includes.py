import nf_core.pipelines.lint

from ..test_lint import TestLint


class TestLintPluginIncludes(TestLint):
    def setUp(self) -> None:
        super().setUp()
        self.new_pipeline = self._make_pipeline_copy()

    def test_default_values_match(self):
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        result = lint_obj.plugin_includes()
        assert len(result["failed"]) == 0
        assert len(result["warned"]) == 0

    def test_wrong_include(self):
        test_path = self.new_pipeline / "test.nf"
        with open(test_path, "w") as of:
            of.write("include { paramsSummary } from 'plugin/nf-validation'\n")
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        result = lint_obj.plugin_includes()
        assert len(result["failed"]) == 1
        assert len(result["warned"]) == 0
