from pathlib import Path

from ruamel.yaml import YAML

import nf_core.pipelines.lint
from nf_core.utils import NFCoreYamlConfig

from ..test_lint import TestLint


class TestLintNfCoreYml(TestLint):
    def setUp(self) -> None:
        super().setUp()
        self.new_pipeline = self._make_pipeline_copy()
        self.nf_core_yml_path = Path(self.new_pipeline) / ".nf-core.yml"
        self.yaml = YAML()
        self.nf_core_yml: NFCoreYamlConfig = self.yaml.load(self.nf_core_yml_path)
        self.lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)

    def test_nfcore_yml_pass(self):
        """Lint test: nfcore_yml - PASS"""
        assert self.lint_obj._load()
        results = self.lint_obj.nfcore_yml()

        assert "Repository type in `.nf-core.yml` is valid" in str(results["passed"])
        assert "nf-core version in `.nf-core.yml` is set to the latest version" in str(results["passed"])
        assert len(results.get("warned", [])) == 0
        assert len(results.get("failed", [])) == 0
        assert len(results.get("ignored", [])) == 0

    def test_nfcore_yml_fail_repo_type(self):
        """Lint test: nfcore_yml - FAIL - repository type not set"""

        self.nf_core_yml["repository_type"] = "foo"
        self.yaml.dump(self.nf_core_yml, self.nf_core_yml_path)
        with self.assertRaises(AssertionError):
            self.lint_obj._load()

    def test_nfcore_yml_fail_nfcore_version(self):
        """Lint test: nfcore_yml - FAIL - nf-core version not set"""

        self.nf_core_yml["nf_core_version"] = "foo"
        self.yaml.dump(self.nf_core_yml, self.nf_core_yml_path)
        assert self.lint_obj._load()
        results = self.lint_obj.nfcore_yml()
        assert "nf-core version in `.nf-core.yml` is not set to the latest version." in str(results["warned"])
        assert len(results.get("failed", [])) == 0
        assert len(results.get("passed", [])) >= 0
        assert len(results.get("ignored", [])) == 0

    def test_nfcore_yml_nested_lint_config(self) -> None:
        """Lint test: nfcore_yml with nested lint config - PASS"""
        valid_yaml = """
        lint:
            files_unchanged:
                - .github/workflows/branch.yml
            # modules_config: False
            modules_config:
                    - fastqc
            # merge_markers: False
            merge_markers:
                    - docs/my_pdf.pdf
            # nextflow_config: False
            nextflow_config:
                - manifest.name
                - config_defaults:
                    - params.annotation_db
                    - params.multiqc_comment_headers
                    - params.custom_table_headers
            multiqc_config:
                - report_section_order
                - report_comment
            files_exist:
                - .github/CONTRIBUTING.md
                - CITATIONS.md
            # template_strings: False
            template_strings:
                    - docs/my_pdf.pdf
        """
        self.nf_core_yml["lint"] = self.yaml.load(valid_yaml)
        self.yaml.dump(self.nf_core_yml, self.nf_core_yml_path)

        assert self.lint_obj._load()
        results = self.lint_obj.nfcore_yml()
        assert len(results.get("failed", [])) == 0
        assert len(results.get("warned", [])) == 0
        assert len(results.get("ignored", [])) == 0

    def test_nfcore_yml_nested_lint_config_bool(self) -> None:
        """Lint test: nfcore_yml with nested lint config - PASS"""
        valid_yaml = """
        lint:
            files_unchanged:
                - .github/workflows/branch.yml
            modules_config: False
            # modules_config:
            #         - fastqc
            merge_markers: False
            # merge_markers:
            #         - docs/my_pdf.pdf
            # nextflow_config: False
            nextflow_config:
                - manifest.name
                - config_defaults:
                    - params.annotation_db
                    - params.multiqc_comment_headers
                    - params.custom_table_headers
            multiqc_config:
                - report_section_order
                - report_comment
            files_exist:
                - .github/CONTRIBUTING.md
                - CITATIONS.md
            template_strings: False
            # template_strings:
            #         - docs/my_pdf.pdf
        """
        self.nf_core_yml["lint"] = self.yaml.load(valid_yaml)
        self.yaml.dump(self.nf_core_yml, self.nf_core_yml_path)

        assert self.lint_obj._load()
        results = self.lint_obj.nfcore_yml()
        assert len(results.get("failed", [])) == 0
        assert len(results.get("warned", [])) == 0
        assert len(results.get("ignored", [])) == 0
