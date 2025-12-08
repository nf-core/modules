from pathlib import Path

import nf_core.pipelines.lint

from ..test_lint import TestLint


class TestLintFilesUnchanged(TestLint):
    def test_files_unchanged_pass(self):
        self.lint_obj._load()
        results = self.lint_obj.files_unchanged()
        assert len(results.get("warned", [])) == 0
        assert len(results.get("failed", [])) == 0
        assert len(results.get("ignored", [])) == 0
        assert not results.get("could_fix", True)

    def test_files_unchanged_fail(self):
        failing_file = Path(".github", "CONTRIBUTING.md")
        new_pipeline = self._make_pipeline_copy()
        with open(Path(new_pipeline, failing_file), "a") as fh:
            fh.write("THIS SHOULD NOT BE HERE")

        lint_obj = nf_core.pipelines.lint.PipelineLint(new_pipeline)
        lint_obj._load()
        results = lint_obj.files_unchanged()
        assert len(results["failed"]) > 0
        assert str(failing_file) in results["failed"][0]
        assert results["could_fix"]
