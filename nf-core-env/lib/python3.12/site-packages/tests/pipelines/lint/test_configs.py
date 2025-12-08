from pathlib import Path

import yaml

import nf_core.pipelines.create
import nf_core.pipelines.lint

from ..test_lint import TestLint


class TestLintConfigs(TestLint):
    def setUp(self) -> None:
        super().setUp()
        self.new_pipeline = self._make_pipeline_copy()

    def test_withname_in_modules_config(self):
        """Tests finding withName in modules.config passes linting."""

        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj._load()
        result = lint_obj.modules_config()
        assert len(result["failed"]) == 0
        assert any(
            ["`FASTQC` found in `conf/modules.config` and Nextflow scripts." in passed for passed in result["passed"]]
        )

    def test_superfluous_withname_in_modules_config_fails(self):
        """Tests finding withName in modules.config fails linting."""

        # Add withName to modules.config
        modules_config = Path(self.new_pipeline) / "conf" / "modules.config"
        with open(modules_config, "a") as f:
            f.write("\nwithName: 'BPIPE' {\n cache = false \n}")
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline, hide_progress=False)
        lint_obj._load()
        result = lint_obj.modules_config()
        assert len(result["failed"]) == 1
        assert result["failed"][0].startswith("`conf/modules.config` contains `withName:BPIPE`")

    def test_ignore_modules_config(self):
        """Tests ignoring the modules.config passes linting."""

        # ignore modules.config in linting
        with open(Path(self.new_pipeline) / ".nf-core.yml") as f:
            content = yaml.safe_load(f)
            old_content = content.copy()
            content["lint"] = {"modules_config": False}
        with open(Path(self.new_pipeline) / ".nf-core.yml", "w") as f:
            yaml.dump(content, f)
        Path(self.new_pipeline, "conf", "modules.config").unlink()
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj._load()
        result = lint_obj.modules_config()
        assert len(result["ignored"]) == 1
        assert result["ignored"][0].startswith("`conf/modules.config` not found, but it is ignored.")
        # cleanup
        with open(Path(self.new_pipeline) / ".nf-core.yml", "w") as f:
            yaml.dump(old_content, f)

    def test_superfluous_withname_in_base_config_fails(self):
        """Tests finding withName in base.config fails linting."""

        # Add withName to base.config
        base_config = Path(self.new_pipeline) / "conf" / "base.config"
        with open(base_config, "a") as f:
            f.write("\nwithName:CUSTOM_DUMPSOFTWAREVERSIONS {\n cache = false \n}")
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj._load()
        result = lint_obj.base_config()
        assert len(result["failed"]) == 1
        assert result["failed"][0].startswith("`conf/base.config` contains `withName:CUSTOM_DUMPSOFTWAREVERSIONS`")

    def test_ignore_base_config(self):
        """Tests ignoring the base.config passes linting."""

        # ignore base.config in linting
        with open(Path(self.new_pipeline) / ".nf-core.yml") as f:
            content = yaml.safe_load(f)
            old_content = content.copy()
            content["lint"] = {"base_config": False}
        with open(Path(self.new_pipeline) / ".nf-core.yml", "w") as f:
            yaml.dump(content, f)
        Path(self.new_pipeline, "conf", "base.config").unlink()
        lint_obj = nf_core.pipelines.lint.PipelineLint(self.new_pipeline)
        lint_obj._load()
        result = lint_obj.base_config()
        assert len(result["ignored"]) == 1
        assert result["ignored"][0].startswith("`conf/base.config` not found, but it is ignored.")
        # cleanup
        with open(Path(self.new_pipeline) / ".nf-core.yml", "w") as f:
            yaml.dump(old_content, f)
