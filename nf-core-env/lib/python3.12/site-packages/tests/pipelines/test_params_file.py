import json
from pathlib import Path

from nf_core.pipelines.params_file import ParamsFileBuilder

from ..test_pipelines import TestPipelines


class TestParamsFileBuilder(TestPipelines):
    """Class for schema tests"""

    def setUp(self):
        """Create a new PipelineSchema object"""
        super().setUp()

        self.template_schema = Path(self.pipeline_dir, "nextflow_schema.json")
        self.params_template_builder = ParamsFileBuilder(self.pipeline_dir)
        self.outfile = Path(self.pipeline_dir, "params-file.yml")

    def test_build_template(self):
        self.params_template_builder.write_params_file(self.outfile)

        assert self.outfile.exists()

        with open(self.outfile) as fh:
            out = fh.read()

        assert f"{self.pipeline_obj.pipeline_prefix}/{self.pipeline_obj.pipeline_name}" in out

    def test_build_template_invalid_schema(self):
        """Build a schema from a template"""
        schema = {}
        with open(self.template_schema) as fh:
            schema = json.load(fh)
            del schema["allOf"]

        with open(self.template_schema, "w") as fh:
            json.dump(schema, fh)

        builder = ParamsFileBuilder(self.template_schema)
        res = builder.write_params_file(self.outfile)

        assert res is False
        assert "Pipeline schema file is invalid" in self.caplog.text

    def test_build_template_file_exists(self):
        """Build a schema from a template"""

        # Creates a new empty file
        self.outfile.touch()

        res = self.params_template_builder.write_params_file(self.outfile)

        assert res is False
        assert f"File '{self.outfile}' exists!" in self.caplog.text

        self.outfile.unlink()

    def test_build_template_content(self):
        """Test that the content of the params file is correct"""
        self.params_template_builder.write_params_file(self.outfile)

        with open(self.outfile) as fh:
            out = fh.read()

        assert f"{self.pipeline_obj.pipeline_prefix}/{self.pipeline_obj.pipeline_name}" in out
        assert "# input: null" in out
