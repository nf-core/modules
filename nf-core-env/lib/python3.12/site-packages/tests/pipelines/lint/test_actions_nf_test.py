from pathlib import Path

import yaml

import nf_core.pipelines.lint

from ..test_lint import TestLint


class TestLintActionsNfTest(TestLint):
    def test_actions_nf_test_pass(self):
        """Lint test: actions_nf_test - PASS"""
        self.lint_obj._load()
        results = self.lint_obj.actions_nf_test()
        assert results["passed"] == [
            "'.github/workflows/nf-test.yml' is triggered on expected events",
            "'.github/workflows/nf-test.yml' checks minimum NF version",
        ]
        assert len(results.get("warned", [])) == 0
        assert len(results.get("failed", [])) == 0
        assert len(results.get("ignored", [])) == 0

    def test_actions_nf_test_fail_wrong_nf(self):
        """Lint test: actions_nf_test - FAIL - wrong minimum version of Nextflow tested"""
        self.lint_obj._load()
        self.lint_obj.minNextflowVersion = "1.2.3"
        results = self.lint_obj.actions_nf_test()
        assert results["failed"] == [
            "Minimum pipeline NF version '1.2.3' is not tested in '.github/workflows/nf-test.yml'"
        ]

    def test_actions_nf_test_fail_wrong_trigger(self):
        """Lint test: actions_actions_nf_test - FAIL - workflow triggered incorrectly, NF ver not checked at all"""

        # Edit .github/workflows/nf-test.yml to mess stuff up!
        new_pipeline = self._make_pipeline_copy()
        with open(Path(new_pipeline, ".github", "workflows", "nf-test.yml")) as fh:
            ci_yml = yaml.safe_load(fh)
        ci_yml[True].pop("pull_request")
        ci_yml["jobs"]["nf-test"]["strategy"]["matrix"] = {"nxf_versionnn": ["foo", ""]}
        with open(Path(new_pipeline, ".github", "workflows", "nf-test.yml"), "w") as fh:
            yaml.dump(ci_yml, fh)

        # Make lint object
        lint_obj = nf_core.pipelines.lint.PipelineLint(new_pipeline)
        lint_obj._load()

        results = lint_obj.actions_nf_test()
        assert results["failed"] == [
            "'.github/workflows/nf-test.yml' is not triggered on expected events",
            "'.github/workflows/nf-test.yml' does not check minimum NF version",
        ]
