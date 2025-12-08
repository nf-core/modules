import nf_core.pipelines.schema


def schema_params(self):
    """Check that the schema describes all flat params in the pipeline.

    The ``nextflow_schema.json`` pipeline schema should describe every flat parameter
    returned from the ``nextflow config`` command (params that are objects or more complex structures are ignored).

    * Failure: If parameters are found in ``nextflow_schema.json`` that are not in ``nextflow_schema.json``
    * Warning: If parameters are found in ``nextflow_schema.json`` that are not in ``nextflow_schema.json``
    """
    passed = []
    warned = []
    failed = []

    # First, get the top-level config options for the pipeline
    self.schema_obj = nf_core.pipelines.schema.PipelineSchema()
    self.schema_obj.get_schema_path(self.wf_path)
    self.schema_obj.get_wf_params()
    self.schema_obj.no_prompts = True
    self.schema_obj.load_lint_schema()

    # Remove any schema params not found in the config
    removed_params = self.schema_obj.remove_schema_notfound_configs()

    # Add schema params found in the config but not the schema
    added_params = self.schema_obj.add_schema_found_configs()

    # Invalid default parameters in nextflow.config
    invalid_config_default_params = self.schema_obj.invalid_nextflow_config_default_parameters

    if len(removed_params) > 0:
        for param in removed_params:
            warned.append(f"Schema param `{param}` not found from nextflow config")

    if len(added_params) > 0:
        for param in added_params:
            failed.append(f"Param `{param}` from `nextflow config` not found in nextflow_schema.json")

    if len(removed_params) == 0 and len(added_params) == 0:
        passed.append("Schema matched params returned from nextflow config")

    if len(invalid_config_default_params) > 0:
        for param, msg in invalid_config_default_params.items():
            failed.append(f"Default value for param `{param}` invalid: {msg}")

    return {"passed": passed, "warned": warned, "failed": failed}
