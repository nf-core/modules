from pathlib import Path

import yaml

import nf_core.pipelines.lint

from ..test_lint import TestLint


class TestLintNfTestContent(TestLint):
    def setUp(self) -> None:
        super().setUp()
        self.new_pipeline = self._make_pipeline_copy()

    def test_nf_test_content_nf_test_files(self):
        """Test failure if nf-test file does not contain outdir parameter."""
        # Create a test file without outdir parameter or versions.yml snapshot
        test_file = Path(self.new_pipeline) / "tests" / "test.nf.test"
        with open(test_file, "w") as f:
            f.write("// No outdir parameter or versions YAML snapshot")
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        result = lint_obj.nf_test_content()
        assert len(result["failed"]) > 0
        assert (
            "'tests/test.nf.test' does not contain `outdir` parameter, it should contain `outdir = \"$outputDir\"`"
            in result["failed"]
        )
        assert "'tests/test.nf.test' does not snapshot a 'versions.yml' file" in result["failed"]

    def test_nf_test_content_nextflow_config_file(self):
        """Test failure if nextflow.config does not contain correct parameters."""
        # Create nextflow.config without required parameters
        config_file = Path(self.new_pipeline) / "tests" / "nextflow.config"
        with open(config_file, "w") as f:
            f.write("// Missing parameters")
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        result = lint_obj.nf_test_content()
        assert len(result["failed"]) > 0
        assert "'tests/nextflow.config' does not contain `modules_testdata_base_path`" in result["failed"]
        assert "'tests/nextflow.config' does not contain `pipelines_testdata_base_path`" in result["failed"]

    def test_nf_test_content_missing_nf_test_config_file(self):
        """Test failure if nf-test.config does not contain required settings."""
        # Create nf-test.config without required settings
        config_file = Path(self.new_pipeline) / "nf-test.config"
        with open(config_file, "w") as f:
            f.write("// Missing required settings")
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        result = lint_obj.nf_test_content()
        assert len(result["failed"]) > 0
        assert "'nf-test.config' does not set a `testsDir`, it should contain `testsDir \".\"`" in result["failed"]
        assert (
            '\'nf-test.config\' does not set a `workDir`, it should contain `workDir System.getenv("NFT_WORKDIR") ?: ".nf-test"`'
            in result["failed"]
        )
        assert (
            "'nf-test.config' does not set a `configFile`, it should contain `configFile \"tests/nextflow.config\"`"
            in result["failed"]
        )

    def test_nf_test_content_ignored(self):
        """Test that nf-test content checks can be ignored via .nf-core.yml."""
        # Create .nf-core.yml to ignore checks
        nf_core_yml = Path(self.new_pipeline) / ".nf-core.yml"
        with open(nf_core_yml) as f:
            yml_content = yaml.safe_load(f)

        yml_content["lint"] = {"nf_test_content": ["tests/default.nf.test", "tests/nextflow.config", "nf-test.config"]}

        with open(nf_core_yml, "w") as f:
            yaml.dump(yml_content, f)

        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj._load()
        result = lint_obj.nf_test_content()
        print(result)
        assert len(result["ignored"]) == 3
        assert "'tests/default.nf.test' checking ignored" in result["ignored"]
        assert "'tests/nextflow.config' checking ignored" in result["ignored"]
        assert "'nf-test.config' checking ignored" in result["ignored"]
