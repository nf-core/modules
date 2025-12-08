from pathlib import Path

import yaml

import nf_core.pipelines.lint

from ..test_lint import TestLint


class TestLintActionsAwsfulltest(TestLint):
    def test_actions_awsfulltest_warn(self):
        """Lint test: actions_awsfulltest - PASS"""
        self.lint_obj._load()
        results = self.lint_obj.actions_awsfulltest()
        assert "`.github/workflows/awsfulltest.yml` is triggered correctly" in results["passed"]
        assert len(results.get("failed", [])) == 0
        assert len(results.get("ignored", [])) == 0

    def test_actions_awsfulltest_pass(self):
        """Lint test: actions_awsfulltest - WARN"""

        # Edit .github/workflows/awsfulltest.yml to use -profile test_full
        new_pipeline = self._make_pipeline_copy()
        with open(Path(new_pipeline, ".github", "workflows", "awsfulltest.yml")) as fh:
            awsfulltest_yml = fh.read()
        awsfulltest_yml = awsfulltest_yml.replace("-profile test ", "-profile test_full ")
        with open(Path(new_pipeline, ".github", "workflows", "awsfulltest.yml"), "w") as fh:
            fh.write(awsfulltest_yml)

        # Make lint object
        lint_obj = nf_core.pipelines.lint.PipelineLint(new_pipeline)
        lint_obj._load()

        results = lint_obj.actions_awsfulltest()
        assert results["passed"] == [
            "`.github/workflows/awsfulltest.yml` is triggered correctly",
            "`.github/workflows/awsfulltest.yml` does not use `-profile test`",
        ]
        assert len(results.get("warned", [])) == 0
        assert len(results.get("failed", [])) == 0
        assert len(results.get("ignored", [])) == 0

    def test_actions_awsfulltest_fail(self):
        """Lint test: actions_awsfulltest - FAIL"""

        # Edit .github/workflows/awsfulltest.yml to use -profile test_full
        new_pipeline = self._make_pipeline_copy()
        with open(Path(new_pipeline, ".github", "workflows", "awsfulltest.yml")) as fh:
            awsfulltest_yml = yaml.safe_load(fh)
        del awsfulltest_yml[True]["pull_request_review"]
        with open(Path(new_pipeline, ".github", "workflows", "awsfulltest.yml"), "w") as fh:
            yaml.dump(awsfulltest_yml, fh)

        # Make lint object
        lint_obj = nf_core.pipelines.lint.PipelineLint(new_pipeline)
        lint_obj._load()

        results = lint_obj.actions_awsfulltest()
        assert results["failed"] == ["`.github/workflows/awsfulltest.yml` is not triggered correctly"]
        assert "`.github/workflows/awsfulltest.yml` does not use `-profile test`" in results["passed"]
        assert len(results.get("ignored", [])) == 0
