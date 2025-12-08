import re
from pathlib import Path

from ruamel.yaml import YAML

import nf_core.pipelines.create.create
import nf_core.pipelines.lint
from nf_core.utils import NFCoreYamlConfig

from ..test_lint import TestLint


class TestLintVersionConsistency(TestLint):
    def setUp(self) -> None:
        super().setUp()
        self.new_pipeline = self._make_pipeline_copy()
        self.nf_core_yml_path = Path(self.new_pipeline) / ".nf-core.yml"
        self.yaml = YAML()
        self.nf_core_yml: NFCoreYamlConfig = self.yaml.load(self.nf_core_yml_path)

    def test_version_consistency_pass(self):
        """Tests that pipeline version consistency test passes with good pipeline version"""
        # Declare the pipeline version in the .nf-core.yml file
        self.nf_core_yml["template"]["version"] = "1.0.0"
        self.yaml.dump(self.nf_core_yml, self.nf_core_yml_path)

        # Set the version in the nextflow.config file
        nf_conf_file = Path(self.new_pipeline, "nextflow.config")
        with open(nf_conf_file) as f:
            content = f.read()
            pass_content = re.sub(r"(?m)^(?P<indent>\s*version\s*=\s*)'[^']+'", rf"\g<indent>'{'1.0.0'}'", content)
        with open(nf_conf_file, "w") as f:
            f.write(pass_content)

        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj.load_pipeline_config()
        lint_obj.nextflow_config()
        # Set the version for the container
        lint_obj.nf_config["process.container"] = "nfcore/pipeline:1.0.0"
        result = lint_obj.version_consistency()
        assert result["passed"] == [
            "Version tags are consistent: manifest.version = 1.0.0, process.container = 1.0.0, nfcore_yml.version = 1.0.0",
        ]
        assert result["failed"] == []

    def test_version_consistency_not_numeric(self):
        """Tests that pipeline version consistency test fails with non-numeric version numbers"""
        self.nf_core_yml["template"]["version"] = "1.0.0dev"
        self.yaml.dump(self.nf_core_yml, self.nf_core_yml_path)

        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj.load_pipeline_config()
        lint_obj.nextflow_config()
        result = lint_obj.version_consistency()
        assert result["passed"] == [
            "Version tags are consistent: manifest.version = 1.0.0dev, nfcore_yml.version = 1.0.0dev"
        ]
        assert result["failed"] == [
            "manifest.version was not numeric: 1.0.0dev!",
            "nfcore_yml.version was not numeric: 1.0.0dev!",
        ]

    def test_version_consistency_not_consistent(self):
        """Tests that pipeline version consistency test fails with inconsistent version numbers"""
        self.nf_core_yml["template"]["version"] = "0.0.0"
        self.yaml.dump(self.nf_core_yml, self.nf_core_yml_path)
        nf_conf_file = Path(self.new_pipeline, "nextflow.config")
        with open(nf_conf_file) as f:
            content = f.read()
            fail_content = re.sub(r"(?m)^(?P<indent>\s*version\s*=\s*)'[^']+'", rf"\g<indent>'{'1.0.0'}'", content)
        with open(nf_conf_file, "w") as f:
            f.write(fail_content)
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj.load_pipeline_config()
        lint_obj.nextflow_config()

        lint_obj.nf_config["process.container"] = "nfcore/pipeline:0.1"
        result = lint_obj.version_consistency()
        assert len(result["passed"]) == 0
        assert result["failed"] == [
            "The versioning is not consistent between container, release tag and config. Found manifest.version = 1.0.0, process.container = 0.1, nfcore_yml.version = 0.0.0",
        ]
