import logging

from nf_core.pipelines.lint.pipeline_todos import pipeline_todos

log = logging.getLogger(__name__)


def module_todos(_, module):
    """
    Look for TODO statements in the module files

    The nf-core module template contains a number of comment lines to help developers
    of new modules know where they need to edit files and add content.
    They typically have the following format:

    .. code-block:: groovy

        // TODO nf-core: Make some kind of change to the workflow here

    ..or in markdown:

    .. code-block:: html

        <!-- TODO nf-core: Add some detail to the docs here -->

    This lint test runs through all files in the module and searches for these lines.
    If any are found they will throw a warning.

    .. tip:: Note that many GUI code editors have plugins to list all instances of *TODO*
              in a given project directory. This is a very quick and convenient way to get
              started on your pipeline!

    """

    # Main module directory
    mod_results = pipeline_todos(None, root_dir=module.component_dir)
    for i, warning in enumerate(mod_results["warned"]):
        module.warned.append(("module_todos", "module_todo", warning, mod_results["file_paths"][i]))
    for i, passed in enumerate(mod_results["passed"]):
        module.passed.append(("module_todos", "module_todo", passed, module.component_dir))
