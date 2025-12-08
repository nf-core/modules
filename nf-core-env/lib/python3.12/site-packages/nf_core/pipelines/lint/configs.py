import logging
import re
from pathlib import Path

from nf_core.pipelines.lint_utils import ignore_file

log = logging.getLogger(__name__)


class LintConfig:
    def __init__(self, wf_path: str, lint_config: dict[str, list[str]]):
        self.wf_path = wf_path
        self.lint_config = lint_config

    def lint_file(self, lint_name: str, file_path: Path) -> dict[str, list[str]]:
        """Lint a file and add the result to the passed or failed list."""

        passed: list[str] = []
        failed: list[str] = []
        ignored: list[str] = []
        ignore_configs: list[str] = []

        fn = Path(self.wf_path, file_path)

        passed, failed, ignored, ignore_configs = ignore_file(lint_name, file_path, Path(self.wf_path))

        error_message = f"`{file_path}` not found"
        # check for partial match in failed or ignored
        if not any(f.startswith(error_message) for f in (failed + ignored)):
            try:
                with open(fn) as fh:
                    config = fh.read()
            except Exception as e:
                return {"failed": [f"Could not parse file: {fn}, {e}"]}

            # find sections with a withName: prefix
            sections = re.findall(r"withName:\s*['\"]?(\w+)['\"]?", config)
            log.debug(f"found sections: {sections}")

            # find all .nf files in the workflow directory
            nf_files = list(Path(self.wf_path).rglob("*.nf"))
            log.debug(f"found nf_files: {[str(f) for f in nf_files]}")

            # check if withName sections are present in config, but not in workflow files
            for section in sections:
                if section not in ignore_configs or section.lower() not in ignore_configs:
                    if not any(section in nf_file.read_text() for nf_file in nf_files):
                        failed.append(
                            f"`{file_path}` contains `withName:{section}`, but the corresponding process is not present in any of the Nextflow scripts."
                        )
                    else:
                        passed.append(f"`{section}` found in `{file_path}` and Nextflow scripts.")
                else:
                    ignored.append(f"``{section}` is ignored.")

        return {"passed": passed, "failed": failed, "ignored": ignored}


def modules_config(self) -> dict[str, list[str]]:
    """Make sure the conf/modules.config file follows the nf-core template, especially removed sections.

    .. note:: You can choose to ignore this lint tests by editing the file called
        ``.nf-core.yml`` in the root of your pipeline and setting the test to false:

        .. code-block:: yaml

            lint:
                modules_config: False

        To disable this test only for specific modules, you can specify a list of module names.

        .. code-block:: yaml

            lint:
                modules_config:
                    - fastqc

    """

    result = LintConfig(self.wf_path, self.lint_config).lint_file("modules_config", Path("conf", "modules.config"))

    return result


def base_config(self) -> dict[str, list[str]]:
    """Make sure the conf/base.config file follows the nf-core template, especially removed sections.

    .. note:: You can choose to ignore this lint tests by editing the file called
        ``.nf-core.yml`` in the root of your pipeline and setting the test to false:

        .. code-block:: yaml

            lint:
                base_config: False

    """

    result = LintConfig(self.wf_path, self.lint_config).lint_file("base_config", Path("conf", "base.config"))

    return result
