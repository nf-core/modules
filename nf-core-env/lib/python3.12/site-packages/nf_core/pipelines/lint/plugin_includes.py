import ast
import glob
import logging
import re

log = logging.getLogger(__name__)


def plugin_includes(self) -> dict[str, list[str]]:
    """Checks the include statements in the all *.nf files for plugin includes

    When nf-schema is used in an nf-core pipeline, the include statements of the plugin
    functions have to use nf-schema instead of nf-validation and vice versa
    """
    config_plugins = [plugin.split("@")[0] for plugin in ast.literal_eval(self.nf_config.get("plugins", "[]"))]
    validation_plugin = "nf-validation" if "nf-validation" in config_plugins else "nf-schema"

    passed: list[str] = []
    warned: list[str] = []
    failed: list[str] = []
    ignored: list[str] = []

    plugin_include_pattern = re.compile(r"^include\s*{[^}]+}\s*from\s*[\"']plugin/([^\"']+)[\"']\s*$", re.MULTILINE)
    workflow_files = [
        file for file in glob.glob(f"{self.wf_path}/**/*.nf", recursive=True) if not file.startswith("./modules/")
    ]
    test_passed = True
    for file in workflow_files:
        with open(file) as of:
            plugin_includes = re.findall(plugin_include_pattern, of.read())
            for include in plugin_includes:
                if include not in ["nf-validation", "nf-schema"]:
                    continue
                if include != validation_plugin:
                    test_passed = False
                    failed.append(
                        f"Found a `{include}` plugin import in `{file[2:]}`, but `{validation_plugin}` was used in `nextflow.config`"
                    )

    if test_passed:
        passed.append("No wrong validation plugin imports have been found")

    return {"passed": passed, "warned": warned, "failed": failed, "ignored": ignored}
