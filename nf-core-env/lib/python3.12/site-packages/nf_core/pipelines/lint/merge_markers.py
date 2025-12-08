import fnmatch
import logging
import os
from pathlib import Path

import nf_core.utils

log = logging.getLogger(__name__)


def merge_markers(self):
    """Check for remaining merge markers.

    This test looks for remaining merge markers in the code, e.g.:
    ``>>>>>>>`` or ``<<<<<<<``

    .. note:: You can choose to ignore this lint tests by editing the file called
        ``.nf-core.yml`` in the root of your pipeline and setting the test to false:

        .. code-block:: yaml

            lint:
                merge_markers: False

        To disable this test only for specific files, you can specify a list of file paths to ignore.
        For example, to ignore a pdf you added to the docs:

        .. code-block:: yaml

            lint:
                merge_markers:
                    - docs/my_pdf.pdf

    """
    passed = []
    failed = []
    ignored = []

    ignored_config = self.lint_config.get("merge_markers", []) if self.lint_config is not None else []

    ignore = [".git"]
    if Path(self.wf_path, ".gitignore").is_file():
        with open(Path(self.wf_path, ".gitignore"), encoding="latin1") as fh:
            for line in fh:
                ignore.append(Path(line.strip().rstrip("/")).name)
    for root, dirs, files in os.walk(self.wf_path, topdown=True):
        # Ignore files
        for i_base in ignore:
            i = str(Path(root, i_base))
            dirs[:] = [d for d in dirs if not fnmatch.fnmatch(str(Path(root, d)), i)]
            files[:] = [f for f in files if not fnmatch.fnmatch(str(Path(root, f)), i)]
        for fname in files:
            # File ignored in config
            if str(Path(root, fname).relative_to(self.wf_path)) in ignored_config:
                ignored.append(f"Ignoring file `{Path(root, fname)}`")
                continue
            # Skip binary files
            if nf_core.utils.is_file_binary(Path(root, fname)):
                continue
            try:
                with open(Path(root, fname), encoding="latin1") as fh:
                    for line in fh:
                        if ">>>>>>>" in line:
                            failed.append(f"Merge marker '>>>>>>>' in `{Path(root, fname)}`: {line[:30]}")
                        if "<<<<<<<" in line:
                            failed.append(f"Merge marker '<<<<<<<' in `{Path(root, fname)}`: {line[:30]}")
            except FileNotFoundError:
                log.debug(f"Could not open file {Path(root, fname)} in merge_markers lint test")
    if len(failed) == 0:
        passed.append("No merge markers found in pipeline files")
    return {"passed": passed, "failed": failed, "ignored": ignored}
