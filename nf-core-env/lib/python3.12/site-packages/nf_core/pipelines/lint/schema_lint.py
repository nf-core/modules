import logging

import nf_core.pipelines.schema


def schema_lint(self):
    """Pipeline schema syntax

    Pipelines should have a ``nextflow_schema.json`` file that describes the different
    pipeline parameters (eg. ``params.something``, ``--something``).

    .. tip:: Reminder: you should generally never need to edit this JSON file by hand.
             The ``nf-core pipelines schema build`` command can create *and edit* the file for you
             to keep it up to date, with a friendly user-interface for customisation.

    The lint test checks the schema for the following:

    * Schema should be a valid JSON file
    * Schema should adhere to `JSONSchema <https://json-schema.org/>`_, Draft 7 or Draft 2020-12.
    * Parameters can be described in two places:

        * As ``properties`` in the top-level schema object
        * As ``properties`` within subschemas listed in a top-level ``definitions`` (draft 7) or ``$defs`` (draft 2020-12) objects

    * The schema must describe at least one parameter
    * There must be no duplicate parameter IDs across the schema and definition subschema
    * All subschema in ``definitions`` or ``$defs`` must be referenced in the top-level ``allOf`` key
    * The top-level ``allOf`` key must not describe any non-existent definitions
    * Default parameters in the schema must be valid
    * Core top-level schema attributes should exist and be set as follows:

        * ``$schema``: ``https://json-schema.org/draft-07/schema`` or ``https://json-schema.org/draft/2020-12/schema``
        * ``$id``: URL to the raw schema file, eg. ``https://raw.githubusercontent.com/YOURPIPELINE/master/nextflow_schema.json``
        * ``title``: ``YOURPIPELINE pipeline parameters``
        * ``description``: The pipeline config ``manifest.description``
    * That the ``input`` property is defined and has a mimetype. A list of common mimetypes can be found `here <https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/MIME_types/Common_types>`_.

    For example, an *extremely* minimal schema could look like this (draft 7):

    .. code-block:: json

       {
         "$schema": "https://json-schema.org/draft-07/schema",
         "$id": "https://raw.githubusercontent.com/YOURPIPELINE/master/nextflow_schema.json",
         "title": "YOURPIPELINE pipeline parameters",
         "description": "This pipeline is for testing",
         "properties": {
           "first_param": { "type": "string" }
         },
         "definitions": {
           "my_first_group": {
             "properties": {
               "second_param": { "type": "string" }
             }
           }
         },
         "allOf": [{"$ref": "#/definitions/my_first_group"}]
       }

    Or this (draft 2020-12):

    .. code-block:: json

       {
         "$schema": "https://json-schema.org/draft/2020-12/schema",
         "$id": "https://raw.githubusercontent.com/YOURPIPELINE/master/nextflow_schema.json",
         "title": "YOURPIPELINE pipeline parameters",
         "description": "This pipeline is for testing",
         "properties": {
           "first_param": { "type": "string" }
         },
         "$defs": {
           "my_first_group": {
             "properties": {
               "second_param": { "type": "string" }
             }
           }
         },
         "allOf": [{"$ref": "#/$defs/my_first_group"}]
       }

    .. tip:: You can check your pipeline schema without having to run the entire pipeline lint
             by running ``nf-core pipelines schema lint`` instead of ``nf-core pipelines lint``
    """
    passed = []
    warned = []
    failed = []

    # Only show error messages from schema
    logging.getLogger("nf_core.pipelines.schema").setLevel(logging.ERROR)

    # Lint the schema
    self.schema_obj = nf_core.pipelines.schema.PipelineSchema()
    self.schema_obj.get_schema_path(self.wf_path)

    try:
        self.schema_obj.load_lint_schema()
        passed.append("Schema lint passed")
    except AssertionError as e:
        failed.append(f"Schema lint failed: {e}")

    # Check the title and description - gives warnings instead of fail
    if self.schema_obj.schema is not None:
        try:
            self.schema_obj.validate_schema_title_description()
            passed.append("Schema title + description lint passed")
        except AssertionError as e:
            warned.append(str(e))

    # Check for mimetype in the 'input' parameter, warn if missing
    if self.schema_obj.schema is not None:
        try:
            has_valid_mimetype = self.schema_obj.check_for_input_mimetype()
            if has_valid_mimetype is not None:
                passed.append(f"Input mimetype lint passed: '{has_valid_mimetype}'")
            else:
                warned.append("Input mimetype is missing or empty")
        except LookupError as e:
            warned.append(str(e))

    return {"passed": passed, "warned": warned, "failed": failed}
