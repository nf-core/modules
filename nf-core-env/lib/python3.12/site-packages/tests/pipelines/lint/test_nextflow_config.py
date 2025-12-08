import os
import re
from pathlib import Path

import yaml

import nf_core.pipelines.create.create
import nf_core.pipelines.lint

from ..test_lint import TestLint


class TestLintNextflowConfig(TestLint):
    def setUp(self) -> None:
        super().setUp()
        self.new_pipeline = self._make_pipeline_copy()

    def test_nextflow_config_example_pass(self):
        """Tests that config variable existence test works with good pipeline example"""
        self.lint_obj.load_pipeline_config()
        result = self.lint_obj.nextflow_config()
        assert len(result["failed"]) == 0
        assert len(result["warned"]) == 0

    def test_default_values_match(self):
        """Test that the default values in nextflow.config match the default values defined in the nextflow_schema.json."""
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj.load_pipeline_config()
        result = lint_obj.nextflow_config()
        assert len(result["failed"]) == 0
        assert len(result["warned"]) == 0
        assert "Config default value correct: params.validate_params" in str(result["passed"])

    def test_nextflow_config_bad_name_fail(self):
        """Tests that config variable existence test fails with bad pipeline name"""
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj.load_pipeline_config()

        lint_obj.nf_config["manifest.name"] = "bad_name"
        result = lint_obj.nextflow_config()
        assert len(result["failed"]) > 0
        assert len(result["warned"]) == 0

    def test_nextflow_config_dev_in_release_mode_failed(self):
        """Tests that config variable existence test fails with dev version in release mode"""
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj.load_pipeline_config()

        lint_obj.release_mode = True
        lint_obj.nf_config["manifest.version"] = "dev_is_bad_name"
        result = lint_obj.nextflow_config()
        assert len(result["failed"]) > 0
        assert len(result["warned"]) == 0

    def test_nextflow_config_missing_test_profile_failed(self):
        """Test failure if config file does not contain `test` profile."""
        # Change the name of the test profile so there is no such profile
        nf_conf_file = os.path.join(self.new_pipeline, "nextflow.config")
        with open(nf_conf_file) as f:
            content = f.read()
            fail_content = re.sub(r"\btest\b", "testfail", content)
        with open(nf_conf_file, "w") as f:
            f.write(fail_content)
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj.load_pipeline_config()
        result = lint_obj.nextflow_config()
        assert len(result["failed"]) > 0
        assert len(result["warned"]) == 0

    def test_default_values_fail(self):
        """Test linting fails if the default values in nextflow.config do not match the ones defined in the nextflow_schema.json."""
        # Change the default value of max_multiqc_email_size in nextflow.config
        nf_conf_file = Path(self.new_pipeline) / "nextflow.config"
        with open(nf_conf_file) as f:
            content = f.read()
            fail_content = re.sub(r"\bmax_multiqc_email_size\s*=\s*'25.MB'", "max_multiqc_email_size = '0'", content)
        with open(nf_conf_file, "w") as f:
            f.write(fail_content)
        # Change the default value of custom_config_version in nextflow_schema.json
        nf_schema_file = Path(self.new_pipeline) / "nextflow_schema.json"
        with open(nf_schema_file) as f:
            content = f.read()
            fail_content = re.sub(r'"default": "master"', '"default": "main"', content)
        with open(nf_schema_file, "w") as f:
            f.write(fail_content)
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj.load_pipeline_config()
        result = lint_obj.nextflow_config()
        assert len(result["failed"]) == 2
        assert (
            "Config default value incorrect: `params.max_multiqc_email_size` is set as `25.MB` in `nextflow_schema.json` but is `0` in `nextflow.config`."
            in result["failed"]
        )
        assert (
            "Config default value incorrect: `params.custom_config_version` is set as `main` in `nextflow_schema.json` but is `master` in `nextflow.config`."
            in result["failed"]
        )

    def test_catch_params_assignment_in_main_nf(self):
        """Test linting fails if main.nf contains an assignment to a parameter from nextflow_schema.json."""
        # Add parameter assignment in main.nf
        main_nf_file = Path(self.new_pipeline) / "main.nf"
        with open(main_nf_file, "a") as f:
            f.write("params.custom_config_base = 'test'")
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj.load_pipeline_config()
        result = lint_obj.nextflow_config()
        assert len(result["failed"]) == 2
        assert (
            result["failed"][1]
            == "Config default value incorrect: `params.custom_config_base` is set as `https://raw.githubusercontent.com/nf-core/configs/master` in `nextflow_schema.json` but is `null` in `nextflow.config`."
        )

    def test_allow_params_reference_in_main_nf(self):
        """Test linting allows for references like `params.aligner == 'bwa'` in main.nf. The test will detect if the bug mentioned in GitHub-issue #2833 reemerges."""
        # Add parameter reference in main.nf
        main_nf_file = Path(self.new_pipeline) / "main.nf"
        with open(main_nf_file, "a") as f:
            f.write("params.custom_config_version == 'main'")
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj.load_pipeline_config()
        result = lint_obj.nextflow_config()
        assert len(result["failed"]) == 0

    def test_default_values_ignored(self):
        """Test ignoring linting of default values."""
        valid_yaml = """
        nextflow_config:
            - manifest.name
            - config_defaults:
                - params.custom_config_version
        """
        # Add custom_config_version to the ignore list
        nf_core_yml_path = Path(self.new_pipeline) / ".nf-core.yml"

        with open(nf_core_yml_path) as f:
            nf_core_yml = yaml.safe_load(f)
            nf_core_yml["lint"] = yaml.safe_load(valid_yaml)
        with open(nf_core_yml_path, "w") as f:
            yaml.dump(nf_core_yml, f)

        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj.load_pipeline_config()
        lint_obj._load_lint_config()
        result = lint_obj.nextflow_config()
        assert len(result["failed"]) == 0
        assert len(result["ignored"]) == 2
        assert "Config default value correct: params.custom_config_version" not in str(result["passed"])
        assert "Config default ignored: params.custom_config_version" in str(result["ignored"])
        assert "Config variable ignored: `manifest.name`" in str(result["ignored"])

    def test_default_values_float(self):
        """Test comparing two float values."""
        # Add a float value `dummy=0.0001` to the nextflow.config below `validate_params`
        nf_conf_file = Path(self.new_pipeline) / "nextflow.config"
        with open(nf_conf_file) as f:
            content = f.read()
            fail_content = re.sub(
                r"validate_params\s*=\s*true",
                "params.validate_params = true\ndummy = 0.000000001",
                content,
            )
        with open(nf_conf_file, "w") as f:
            f.write(fail_content)
        # Add a float value `dummy` to the nextflow_schema.json
        nf_schema_file = Path(self.new_pipeline) / "nextflow_schema.json"
        with open(nf_schema_file) as f:
            content = f.read()
            fail_content = re.sub(
                r'"validate_params": {',
                '    "dummy": {"type": "number","default":0.000000001},\n"validate_params": {',
                content,
            )
        with open(nf_schema_file, "w") as f:
            f.write(fail_content)

        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj.load_pipeline_config()
        result = lint_obj.nextflow_config()
        assert len(result["failed"]) == 0
        assert len(result["warned"]) == 0
        assert "Config default value correct: params.dummy" in str(result["passed"])

    def test_default_values_float_fail(self):
        """Test comparing two float values."""
        # Add a float value `dummy=0.0001` to the nextflow.config below `validate_params`
        nf_conf_file = Path(self.new_pipeline) / "nextflow.config"
        with open(nf_conf_file) as f:
            content = f.read()
            fail_content = re.sub(
                r"validate_params\s*=\s*true",
                "params.validate_params = true\ndummy = 0.000000001",
                content,
            )
        with open(nf_conf_file, "w") as f:
            f.write(fail_content)
        # Add a float value `dummy` to the nextflow_schema.json
        nf_schema_file = Path(self.new_pipeline) / "nextflow_schema.json"
        with open(nf_schema_file) as f:
            content = f.read()
            fail_content = re.sub(
                r'"validate_params": {',
                '    "dummy": {"type": "float","default":0.000001},\n"validate_params": {',
                content,
            )
        with open(nf_schema_file, "w") as f:
            f.write(fail_content)

        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj.load_pipeline_config()
        result = lint_obj.nextflow_config()

        assert len(result["failed"]) == 1
        assert len(result["warned"]) == 0
        assert "Config default value incorrect: `params.dummy" in str(result["failed"])
