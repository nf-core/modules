import filecmp
import logging
import os
import re
import shutil
import tempfile
from pathlib import Path

import yaml

import nf_core.pipelines.create.create

log = logging.getLogger(__name__)


def files_unchanged(self) -> dict[str, list[str] | bool]:
    """Checks that certain pipeline files are not modified from template output.

    Iterates through the pipeline's directory content and compares specified files
    against output from the template using the pipeline's metadata. File content
    should not be modified / missing.

    Files that must be unchanged::

        .gitattributes
        .prettierrc.yml
        .github/.dockstore.yml
        .github/CONTRIBUTING.md
        .github/ISSUE_TEMPLATE/bug_report.yml
        .github/ISSUE_TEMPLATE/config.yml
        .github/ISSUE_TEMPLATE/feature_request.yml
        .github/PULL_REQUEST_TEMPLATE.md
        .github/workflows/branch.yml
        .github/workflows/linting_comment.yml
        .github/workflows/linting.yml
        assets/email_template.html
        assets/email_template.txt
        assets/nf-core-PIPELINE_logo_light.png
        assets/sendmail_template.txt
        CODE_OF_CONDUCT.md
        docs/images/nf-core-PIPELINE_logo_light.png
        docs/images/nf-core-PIPELINE_logo_dark.png
        docs/README.md'
        ['LICENSE', 'LICENSE.md', 'LICENCE', 'LICENCE.md'], # NB: British / American spelling

    Files that can have additional content but must include the template contents::

        .gitignore
        .prettierignore


    .. tip:: You can configure the ``nf-core pipelines lint`` tests to ignore any of these checks by setting
             the ``files_unchanged`` key as follows in your ``.nf-core.yml`` config file. For example:

             .. code-block:: yaml

                lint:
                    files_unchanged:
                        - .github/workflows/branch.yml

    """

    passed: list[str] = []
    failed: list[str] = []
    warned: list[str] = []
    ignored: list[str] = []
    fixed: list[str] = []
    could_fix: bool = False

    # Check that we have the minimum required config
    required_pipeline_config = {
        "manifest.name",
        "manifest.description",
        "manifest.contributors",
    }
    missing_pipeline_config = required_pipeline_config.difference(self.nf_config)
    if missing_pipeline_config:
        return {"ignored": [f"Required pipeline config not found - {missing_pipeline_config}"]}
    try:
        prefix, short_name = self.nf_config["manifest.name"].strip("\"'").split("/")
    except ValueError:
        log.warning(
            "Expected manifest.name to be in the format '<repo>/<pipeline>'. Will assume it is <pipeline> and default to repo 'nf-core'"
        )
        short_name = self.nf_config["manifest.name"].strip("\"'")
        prefix = "nf-core"

    # NB: Should all be files, not directories
    # List of lists. Passes if any of the files in the sublist are found.
    files_exact = [
        [Path(".gitattributes")],
        [Path(".prettierrc.yml")],
        [Path("CODE_OF_CONDUCT.md")],
        [Path("LICENSE"), Path("LICENSE.md"), Path("LICENCE"), Path("LICENCE.md")],  # NB: British / American spelling
        [Path(".github", ".dockstore.yml")],
        [Path(".github", "CONTRIBUTING.md")],
        [Path(".github", "ISSUE_TEMPLATE", "bug_report.yml")],
        [Path(".github", "ISSUE_TEMPLATE", "config.yml")],
        [Path(".github", "ISSUE_TEMPLATE", "feature_request.yml")],
        [Path(".github", "PULL_REQUEST_TEMPLATE.md")],
        [Path(".github", "workflows", "branch.yml")],
        [Path(".github", "workflows", "linting_comment.yml")],
        [Path(".github", "workflows", "linting.yml")],
        [Path("assets", "email_template.html")],
        [Path("assets", "email_template.txt")],
        [Path("assets", "sendmail_template.txt")],
        [Path("assets", f"nf-core-{short_name}_logo_light.png")],
        [Path("docs", "images", f"nf-core-{short_name}_logo_light.png")],
        [Path("docs", "images", f"nf-core-{short_name}_logo_dark.png")],
        [Path("docs", "README.md")],
    ]
    files_partial = [
        [Path(".gitignore"), Path(".prettierignore")],
    ]

    # Only show error messages from pipeline creation
    logging.getLogger("nf_core.pipelines.create").setLevel(logging.ERROR)

    # Generate a new pipeline with nf-core create that we can compare to
    tmp_dir = Path(tempfile.TemporaryDirectory().name)
    tmp_dir.mkdir(parents=True)

    # Create a template.yaml file for the pipeline creation
    if "manifest.author" in self.nf_config:
        names = self.nf_config["manifest.author"].strip("\"'")
    if "manifest.contributors" in self.nf_config:
        contributors = self.nf_config["manifest.contributors"]
        names = ", ".join(re.findall(r"name:'([^']+)'", contributors))
    template_yaml = {
        "name": short_name,
        "description": self.nf_config["manifest.description"].strip("\"'"),
        "author": names,
        "org": prefix,
    }

    template_yaml_path = Path(tmp_dir, "template.yaml")

    with open(template_yaml_path, "w") as fh:
        yaml.dump(template_yaml, fh, default_flow_style=False)

    test_pipeline_dir = Path(tmp_dir, f"{prefix}-{short_name}")
    create_obj = nf_core.pipelines.create.create.PipelineCreate(
        None, None, None, no_git=True, outdir=test_pipeline_dir, template_config=template_yaml_path
    )
    create_obj.init_pipeline()

    # Helper functions for file paths
    def _pf(file_path: str | Path) -> Path:
        """Helper function - get file path for pipeline file"""
        return Path(self.wf_path, file_path)

    def _tf(file_path: str | Path) -> Path:
        """Helper function - get file path for template file"""
        return Path(test_pipeline_dir, file_path)

    ignore_files = self.lint_config.get("files_unchanged", []) if self.lint_config is not None else []

    # Files that must be completely unchanged from template
    for files in files_exact:
        # Ignore if file specified in linting config
        if any([str(f) in ignore_files for f in files]):
            ignored.append(f"File ignored due to lint config: {self._wrap_quotes(files)}")

        # Ignore if we can't find the file
        elif not any([_pf(f).is_file() for f in files]):
            ignored.append(f"File does not exist: {self._wrap_quotes(files)}")

        # Check that the file has an identical match
        else:
            for f in files:
                try:
                    if filecmp.cmp(_pf(f), _tf(f), shallow=True):
                        passed.append(f"`{f}` matches the template")
                    else:
                        if f.name.endswith(".png") and int(os.stat(_pf(f)).st_size / 500) == int(
                            os.stat(_tf(f)).st_size / 500
                        ):
                            # almost the same file, good enough for the logo
                            log.debug(f"Files are almost the same. Will pass: {f}")
                            passed.append(f"`{f}` matches the template")
                        elif "files_unchanged" in self.fix:
                            # Try to fix the problem by overwriting the pipeline file
                            shutil.copy(_tf(f), _pf(f))
                            passed.append(f"`{f}` matches the template")
                            fixed.append(f"`{f}` overwritten with template file")
                        elif f.name in ["LICENSE", "LICENSE.md", "LICENCE", "LICENCE.md"]:
                            # Report LICENSE as a warning since we are not using the manifest.author names
                            # TODO: Lint the content of the LICENSE file except the line containing author names
                            # to allow for people to opt-in listing author/maintainer names instead of using the "nf-core community"
                            warned.append(f"`{f}` does not match the template")
                            could_fix = True
                        else:
                            failed.append(f"`{f}` does not match the template")
                            could_fix = True
                except FileNotFoundError:
                    pass

    # Files that can be added to, but that must contain the template contents
    for files in files_partial:
        # Ignore if file specified in linting config
        if any([str(f) in ignore_files for f in files]):
            ignored.append(f"File ignored due to lint config: {self._wrap_quotes(files)}")

        # Ignore if we can't find the file
        elif not any([_pf(f).is_file() for f in files]):
            ignored.append(f"File does not exist: {self._wrap_quotes(files)}")

        # Check that the file contains the template file contents
        else:
            for f in files:
                try:
                    with open(_pf(f)) as fh:
                        pipeline_file = fh.read()
                    with open(_tf(f)) as fh:
                        template_file = fh.read()
                    if template_file in pipeline_file:
                        passed.append(f"`{f}` matches the template")
                    else:
                        if "files_unchanged" in self.fix:
                            # Try to fix the problem by overwriting the pipeline file
                            with open(_tf(f)) as fh:
                                template_file = fh.read()
                            with open(_pf(f), "w") as fh:
                                fh.write(template_file)
                            passed.append(f"`{f}` matches the template")
                            fixed.append(f"`{f}` overwritten with template file")
                        else:
                            failed.append(f"`{f}` does not match the template")
                            could_fix = True
                except FileNotFoundError:
                    pass

    # cleaning up temporary dir
    shutil.rmtree(tmp_dir)

    return {
        "passed": passed,
        "failed": failed,
        "warned": warned,
        "ignored": ignored,
        "fixed": fixed,
        "could_fix": could_fix,
    }
