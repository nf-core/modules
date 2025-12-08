import logging

from nf_core.pipelines.lint.pipeline_todos import pipeline_todos

log = logging.getLogger(__name__)


def subworkflow_todos(_, subworkflow):
    """
    Look for TODO statements in the subworkflow files

    The nf-core subworkflow template contains a number of comment lines to help developers
    of new subworkflow know where they need to edit files and add content.
    They typically have the following format:

    .. code-block:: groovy

        // TODO nf-core: Make some kind of change to the workflow here

    ..or in markdown:

    .. code-block:: html

        <!-- TODO nf-core: Add some detail to the docs here -->

    This lint test runs through all files in the subworkflows and searches for these lines.
    If any are found they will throw a warning.

    .. tip:: Note that many GUI code editors have plugins to list all instances of *TODO*
              in a given project directory. This is a very quick and convenient way to get
              started on your pipeline!

    """

    # Main subworkflow directory
    swf_results = pipeline_todos(None, root_dir=subworkflow.component_dir)
    for i, warning in enumerate(swf_results["warned"]):
        subworkflow.warned.append(("subworkflow_todos", "subworkflow_todo", warning, swf_results["file_paths"][i]))
    for i, passed in enumerate(swf_results["passed"]):
        subworkflow.passed.append(("subworkflow_todos", "subworkflow_todo", passed, subworkflow.component_dir))
