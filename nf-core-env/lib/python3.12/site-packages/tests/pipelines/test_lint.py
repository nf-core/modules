"""Some tests covering the linting code."""

import json
from pathlib import Path

import yaml

import nf_core.pipelines.create.create
import nf_core.pipelines.lint

from ..test_pipelines import TestPipelines
from ..utils import with_temporary_folder


class TestLint(TestPipelines):
    """Class for lint tests"""

    def setUp(self) -> None:
        super().setUp()
        self.lint_obj = nf_core.pipelines.lint.PipelineLint(self.pipeline_dir)


##########################
# CORE lint.py FUNCTIONS #
##########################
class TestPipelinesLint(TestLint):
    def test_run_linting_function(self):
        """Run the run_linting() function in lint.py

        We don't really check any of this code as it's just a series of function calls
        and we're testing each of those individually. This is mostly to check for syntax errors."""
        nf_core.pipelines.lint.run_linting(self.pipeline_dir, False)

    def test_init_pipeline_lint(self):
        """Simply create a PipelineLint object.

        This checks that all of the lint test imports are working properly,
        we also check that the git sha was found and that the release flag works properly
        """
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.pipeline_dir, True)

        # Tests that extra test is added for release mode
        assert "version_consistency" in lint_obj.lint_tests
        assert lint_obj.git_sha
        # Tests that parent nf_core.utils.Pipeline class __init__() is working to find git hash
        assert len(lint_obj.git_sha) > 0

    def test_load_lint_config_not_found(self):
        """Try to load a linting config file that doesn't exist"""
        assert self.lint_obj._load_lint_config()
        assert self.lint_obj.lint_config is not None
        assert self.lint_obj.lint_config.model_dump(exclude_none=True) == {}

    def test_load_lint_config_ignore_all_tests(self):
        """Try to load a linting config file that ignores all tests"""

        # Make a copy of the test pipeline and create a lint object
        new_pipeline = self._make_pipeline_copy()
        lint_obj = nf_core.pipelines.lint.PipelineLint(new_pipeline)

        # Make a config file listing all test names
        config_dict = {"repository_type": "pipeline", "lint": {test_name: False for test_name in lint_obj.lint_tests}}
        with open(Path(new_pipeline, ".nf-core.yml"), "w") as fh:
            yaml.dump(config_dict, fh)

        # Load the new lint config file and check
        lint_obj._load_lint_config()
        assert lint_obj.lint_config is not None
        assert sorted(list(lint_obj.lint_config.model_dump(exclude_none=True))) == sorted(lint_obj.lint_tests)

        # Try running linting and make sure that all tests are ignored
        lint_obj._lint_pipeline()
        assert len(lint_obj.passed) == 0
        assert len(lint_obj.warned) == 0
        assert len(lint_obj.failed) == 0
        assert len(lint_obj.ignored) == len(lint_obj.lint_tests)

    @with_temporary_folder
    def test_json_output(self, tmp_dir):
        """
        Test creation of a JSON file with lint results

        Expected JSON output:
        {
            "nf_core_tools_version": "1.10.dev0",
            "date_run": "2020-06-05 10:56:42",
            "tests_pass": [
                [ 1, "This test passed"],
                [ 2, "This test also passed"]
            ],
            "tests_warned": [
                [ 2, "This test gave a warning"]
            ],
            "tests_failed": [],
            "num_tests_pass": 2,
            "num_tests_warned": 1,
            "num_tests_failed": 0,
            "has_tests_pass": true,
            "has_tests_warned": true,
            "has_tests_failed": false
        }
        """
        self.lint_obj.passed.append(("test_one", "This test passed"))
        self.lint_obj.passed.append(("test_two", "This test also passed"))
        self.lint_obj.warned.append(("test_three", "This test gave a warning"))

        # Make a temp dir for the JSON output
        json_fn = Path(tmp_dir, "lint_results.json")
        self.lint_obj._save_json_results(json_fn)

        # Load created JSON file and check its contents
        with open(json_fn) as fh:
            try:
                saved_json = json.load(fh)
            except json.JSONDecodeError as e:
                raise UserWarning(f"Unable to load JSON file '{json_fn}' due to error {e}")
        assert saved_json["num_tests_pass"] > 0
        assert saved_json["num_tests_warned"] > 0
        assert saved_json["num_tests_ignored"] == 0
        assert saved_json["num_tests_failed"] == 0
        assert saved_json["has_tests_pass"]
        assert saved_json["has_tests_warned"]
        assert not saved_json["has_tests_ignored"]
        assert not saved_json["has_tests_failed"]

    def test_wrap_quotes(self):
        md = self.lint_obj._wrap_quotes(["one", "two", "three"])
        assert md == "`one` or `two` or `three`"

    def test_sphinx_md_files(self):
        """Check that we have .md files for all lint module code,
        and that there are no unexpected files (eg. deleted lint tests)
        To auto-generate the doc files for linting, you can run: `python docs/api/make_lint_md.py`"""

        docs_basedir = Path(Path(__file__).parent.parent.parent, "docs", "api", "_src", "pipeline_lint_tests")

        # Get list of existing .md files
        existing_docs = []
        existing_docs = [
            str(Path(docs_basedir, fn))
            for fn in Path(docs_basedir).iterdir()
            if fn.match("*.md") and not fn.match("index.md")
        ]

        # Check .md files against each test name
        lint_obj = nf_core.pipelines.lint.PipelineLint("", True)
        for test_name in lint_obj.lint_tests:
            fn = Path(docs_basedir, f"{test_name}.md")
            assert fn.exists(), f"Could not find lint docs .md file: {fn}"
            existing_docs.remove(str(fn))

        # Check that we have no remaining .md files that we didn't expect
        assert len(existing_docs) == 0, f"Unexpected lint docs .md files found: {', '.join(existing_docs)}"
