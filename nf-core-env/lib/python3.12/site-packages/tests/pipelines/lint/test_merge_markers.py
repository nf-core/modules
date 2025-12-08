import os

import nf_core.pipelines.lint

from ..test_lint import TestLint


class TestLintMergeMarkers(TestLint):
    def test_merge_markers_found(self):
        """Missing 'jobs' field should result in failure"""
        new_pipeline = self._make_pipeline_copy()

        with open(os.path.join(new_pipeline, "main.nf")) as fh:
            main_nf_content = fh.read()
        main_nf_content = ">>>>>>>\n" + main_nf_content
        with open(os.path.join(new_pipeline, "main.nf"), "w") as fh:
            fh.write(main_nf_content)

        lint_obj = nf_core.pipelines.lint.PipelineLint(new_pipeline)
        lint_obj._load()

        results = lint_obj.merge_markers()
        assert len(results["failed"]) > 0
        assert len(results["passed"]) == 0
        assert "Merge marker '>>>>>>>' in " in results["failed"][0]
