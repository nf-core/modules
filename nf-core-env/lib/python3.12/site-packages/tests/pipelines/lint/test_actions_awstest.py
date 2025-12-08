from pathlib import Path

import yaml

import nf_core.pipelines.lint

from ..test_lint import TestLint


class TestLintActionsAws(TestLint):
    def test_actions_awstest_pass(self):
        """Lint test: actions_awstest - PASS"""
        self.lint_obj._load()
        results = self.lint_obj.actions_awstest()
        assert results["passed"] == ["'.github/workflows/awstest.yml' is triggered correctly"]
        assert len(results.get("warned", [])) == 0
        assert len(results.get("failed", [])) == 0
        assert len(results.get("ignored", [])) == 0

    def test_actions_awstest_fail(self):
        """Lint test: actions_awsfulltest - FAIL"""

        # Edit .github/workflows/awsfulltest.yml to use -profile test_full
        new_pipeline = self._make_pipeline_copy()
        with open(Path(new_pipeline, ".github", "workflows", "awstest.yml")) as fh:
            awstest_yml = yaml.safe_load(fh)
        awstest_yml[True]["push"] = ["main"]
        with open(Path(new_pipeline, ".github", "workflows", "awstest.yml"), "w") as fh:
            yaml.dump(awstest_yml, fh)

        # Make lint object
        lint_obj = nf_core.pipelines.lint.PipelineLint(new_pipeline)
        lint_obj._load()

        results = lint_obj.actions_awstest()
        assert results["failed"] == ["'.github/workflows/awstest.yml' is not triggered correctly"]
        assert len(results.get("warned", [])) == 0
        assert len(results.get("passed", [])) == 0
        assert len(results.get("ignored", [])) == 0
