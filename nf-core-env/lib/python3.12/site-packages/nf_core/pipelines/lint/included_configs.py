from pathlib import Path


def included_configs(self):
    """Check that the pipeline nextflow.config includes the pipeline custom configs.

    If the include line is uncommented, the test passes.
    If the include line is commented, the test fails.
    If the include line is missing, the test warns.

    Can be skipped by adding the following to the .nf-core.yml file:
    lint:
        included_configs: False
    """
    passed = []
    failed = []
    warned = []

    config_file = Path(self.wf_path / "nextflow.config")

    with open(config_file) as fh:
        config = fh.read()
        if (
            f"// includeConfig params.custom_config_base && (!System.getenv('NXF_OFFLINE') || !params.custom_config_base.startsWith('http')) ? \"${{params.custom_config_base}}/pipeline/{self.pipeline_name}.config\""
            in config
        ):
            failed.append("Pipeline config does not include custom configs. Please uncomment the includeConfig line.")
        elif (
            f"includeConfig params.custom_config_base && (!System.getenv('NXF_OFFLINE') || !params.custom_config_base.startsWith('http')) ? \"${{params.custom_config_base}}/pipeline/{self.pipeline_name}.config\""
            in config
        ):
            passed.append("Pipeline config includes custom configs.")
        else:
            warned.append("Pipeline config does not include custom configs. Please add the includeConfig line.")

    return {"passed": passed, "failed": failed, "warned": warned}
