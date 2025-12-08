import logging

from nf_core.pipelines.lint.pipeline_if_empty_null import pipeline_if_empty_null

log = logging.getLogger(__name__)


def subworkflow_if_empty_null(_, subworkflow):
    """Check for ifEmpty(null)

    There are two general cases for workflows to use the channel operator `ifEmpty`:
        1. `ifEmpty( [ ] )` to ensure a process executes, for example when an input file is optional (although this can be replaced by `toList()`).
        2. When a channel should not be empty and throws an error `ifEmpty { error ... }`, e.g. reading from an empty samplesheet.

    There are multiple examples of workflows that inject null objects into channels using `ifEmpty(null)`, which can cause unhandled null pointer exceptions.
    This lint test throws warnings for those instances.
    """

    # Main subworkflow directory
    swf_results = pipeline_if_empty_null(None, root_dir=subworkflow.component_dir)
    for i, warning in enumerate(swf_results["warned"]):
        subworkflow.warned.append(
            ("subworkflow_if_empty_null", "subworkflow_if_empty_null", warning, swf_results["file_paths"][i])
        )
    for i, passed in enumerate(swf_results["passed"]):
        subworkflow.passed.append(
            ("subworkflow_if_empty_null", "subworkflow_if_empty_null", passed, subworkflow.component_dir)
        )
