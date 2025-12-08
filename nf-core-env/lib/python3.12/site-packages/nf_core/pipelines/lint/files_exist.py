import logging
from pathlib import Path

log = logging.getLogger(__name__)


def files_exist(self) -> dict[str, list[str]]:
    """Checks a given pipeline directory for required files.

    Iterates through the pipeline's directory content and checks that specified
    files are either present or absent, as required.

    .. note::
        This test raises an ``AssertionError`` if neither ``nextflow.config`` or ``main.nf`` are found.
        If these files are not found then this cannot be a Nextflow pipeline and something has gone badly wrong.
        All lint tests are stopped immediately with a critical error message.

    Files that *must* be present:

    .. code-block:: bash

        .gitattributes
        .gitignore
        .nf-core.yml
        .prettierignore
        .prettierrc.yml
        .github/.dockstore.yml
        .github/CONTRIBUTING.md
        .github/ISSUE_TEMPLATE/bug_report.yml
        .github/ISSUE_TEMPLATE/config.yml
        .github/ISSUE_TEMPLATE/feature_request.yml
        .github/PULL_REQUEST_TEMPLATE.md
        .github/workflows/branch.yml
        .github/workflows/nf-test.yml
        .github/actions/get-shards/action.yml
        .github/actions/nf-test/action.yml
        .github/workflows/linting_comment.yml
        .github/workflows/linting.yml
        [LICENSE, LICENSE.md, LICENCE, LICENCE.md]  # NB: British / American spelling
        assets/email_template.html
        assets/email_template.txt
        assets/nf-core-PIPELINE_logo_light.png
        assets/sendmail_template.txt
        conf/modules.config
        conf/test.config
        conf/test_full.config
        CHANGELOG.md
        CITATIONS.md
        CODE_OF_CONDUCT.md
        docs/images/nf-core-PIPELINE_logo_light.png
        docs/images/nf-core-PIPELINE_logo_dark.png
        docs/output.md
        docs/README.md
        docs/usage.md
        nextflow_schema.json
        nextflow.config
        nf-test.config
        README.md
        tests/default.nf.test

    Files that *should* be present:

    .. code-block:: bash

        main.nf
        assets/multiqc_config.yml
        conf/base.config
        conf/igenomes.config
        .github/workflows/awstest.yml
        .github/workflows/awsfulltest.yml
        ro-crate-metadata.json

    Files that *must not* be present, due to being renamed or removed in the template:

    .. code-block:: bash

        .github/ISSUE_TEMPLATE/bug_report.md
        .github/ISSUE_TEMPLATE/feature_request.md
        .github/workflows/push_dockerhub.yml
        .markdownlint.yml
        .nf-core.yaml  # NB: Should be yml, not yaml
        .yamllint.yml
        bin/markdown_to_html.r
        conf/aws.config
        docs/images/nf-core-PIPELINE_logo.png
        lib/Checks.groovy
        lib/Completion.groovy
        lib/NfcoreTemplate.groovy
        lib/Utils.groovy
        lib/Workflow.groovy
        lib/WorkflowMain.groovy
        lib/WorkflowPIPELINE.groovy
        lib/nfcore_external_java_deps.jar
        parameters.settings.json
        pipeline_template.yml # saving information in .nf-core.yml
        Singularity


    Files that *should not* be present:

    .. code-block:: bash

        .travis.yml

    .. tip:: You can configure the ``nf-core pipelines lint`` tests to ignore any of these checks by setting
            the ``files_exist`` key as follows in your ``.nf-core.yml`` config file. For example:

            .. code-block:: yaml

            lint:
                files_exist:
                    - assets/multiqc_config.yml
    """

    passed = []
    warned = []
    failed = []
    ignored = []

    # NB: Should all be files, not directories
    # List of lists. Passes if any of the files in the sublist are found.
    #: test autodoc
    try:
        _, short_name = self.nf_config["manifest.name"].strip("\"'").split("/")
    except ValueError:
        log.warning("Expected manifest.name to be in the format '<repo>/<pipeline>'. Will assume it is '<pipeline>'.")
        short_name = self.nf_config["manifest.name"].strip("\"'").split("/")

    files_fail = [
        [Path(".gitattributes")],
        [Path(".gitignore")],
        [Path(".nf-core.yml")],
        [Path(".prettierignore")],
        [Path(".prettierrc.yml")],
        [Path("CHANGELOG.md")],
        [Path("CITATIONS.md")],
        [Path("CODE_OF_CONDUCT.md")],
        [Path("LICENSE"), Path("LICENSE.md"), Path("LICENCE"), Path("LICENCE.md")],  # NB: British / American spelling
        [Path("nextflow_schema.json")],
        [Path("nextflow.config")],
        [Path("README.md")],
        [Path(".github", ".dockstore.yml")],
        [Path(".github", "CONTRIBUTING.md")],
        [Path(".github", "ISSUE_TEMPLATE", "bug_report.yml")],
        [Path(".github", "ISSUE_TEMPLATE", "config.yml")],
        [Path(".github", "ISSUE_TEMPLATE", "feature_request.yml")],
        [Path(".github", "PULL_REQUEST_TEMPLATE.md")],
        [Path(".github", "workflows", "branch.yml")],
        [Path(".github", "workflows", "nf-test.yml")],
        [Path(".github", "actions", "get-shards", "action.yml")],
        [Path(".github", "actions", "nf-test", "action.yml")],
        [Path(".github", "workflows", "linting_comment.yml")],
        [Path(".github", "workflows", "linting.yml")],
        [Path("assets", "email_template.html")],
        [Path("assets", "email_template.txt")],
        [Path("assets", "sendmail_template.txt")],
        [Path("assets", f"nf-core-{short_name}_logo_light.png")],
        [Path("conf", "modules.config")],
        [Path("conf", "test.config")],
        [Path("conf", "test_full.config")],
        [Path("docs", "images", f"nf-core-{short_name}_logo_light.png")],
        [Path("docs", "images", f"nf-core-{short_name}_logo_dark.png")],
        [Path("docs", "output.md")],
        [Path("docs", "README.md")],
        [Path("docs", "README.md")],
        [Path("docs", "usage.md")],
        [Path("nf-test.config")],
        [Path("tests", "default.nf.test")],
    ]

    files_warn = [
        [Path("main.nf")],
        [Path("assets", "multiqc_config.yml")],
        [Path("conf", "base.config")],
        [Path("conf", "igenomes.config")],
        [Path("conf", "igenomes_ignored.config")],
        [Path(".github", "workflows", "awstest.yml")],
        [Path(".github", "workflows", "awsfulltest.yml")],
        [Path("modules.json")],
        [Path("ro-crate-metadata.json")],
    ]

    # List of strings. Fails / warns if any of the strings exist.
    files_fail_ifexists = [
        Path(".github", "ISSUE_TEMPLATE", "bug_report.md"),
        Path(".github", "ISSUE_TEMPLATE", "feature_request.md"),
        Path(".github", "workflows", "push_dockerhub.yml"),
        Path(".markdownlint.yml"),
        Path(".nf-core.yaml"),  # yml not yaml
        Path(".yamllint.yml"),
        Path("bin", "markdown_to_html.r"),
        Path("conf", "aws.config"),
        Path("docs", "images", f"nf-core-{short_name}_logo.png"),
        Path("lib", "Checks.groovy"),
        Path("lib", "Completion.groovy"),
        Path("lib", "NfcoreTemplate.groovy"),
        Path("lib", "Utils.groovy"),
        Path("lib", "Workflow.groovy"),
        Path("lib", "WorkflowMain.groovy"),
        Path("lib", f"Workflow{short_name[0].upper()}{short_name[1:]}.groovy"),
        Path("parameters.settings.json"),
        Path("pipeline_template.yml"),  # saving information in .nf-core.yml
        Path("Singularity"),
        Path("lib", "nfcore_external_java_deps.jar"),
    ]
    files_warn_ifexists = [Path(".travis.yml")]

    files_hint = [
        [
            ["ro-crate-metadata.json"],
            ". Run `nf-core rocrate` to generate this file. Read more about RO-Crates in the [nf-core/tools docs](https://nf-co.re/tools#create-a-ro-crate-metadata-file).",
        ],
    ]
    # Remove files that should be ignored according to the linting config
    ignore_files = self.lint_config.get("files_exist", []) if self.lint_config is not None else []

    def pf(file_path: str | Path) -> Path:
        return Path(self.wf_path, file_path)

    # First - critical files. Check that this is actually a Nextflow pipeline
    if not pf("nextflow.config").is_file() and not pf("main.nf").is_file():
        failed.append("File not found: nextflow.config or main.nf")
        raise AssertionError("Neither nextflow.config or main.nf found! Is this a Nextflow pipeline?")

    # Files that cause an error if they don't exist
    for files in files_fail:
        if any([str(f) in ignore_files for f in files]):
            continue
        if any([pf(f).is_file() for f in files]):
            passed.append(f"File found: {self._wrap_quotes(files)}")
        else:
            failed.append(f"File not found: {self._wrap_quotes(files)}")

    # Files that cause a warning if they don't exist
    for files in files_warn:
        if any([str(f) in ignore_files for f in files]):
            continue
        if any([pf(f).is_file() for f in files]):
            passed.append(f"File found: {self._wrap_quotes(files)}")
        else:
            hint = ""
            for file_hint in files_hint:
                if file_hint[0] == files:
                    hint = str(file_hint[1])
            warned.append(f"File not found: {self._wrap_quotes(files)}{hint}")

    # Files that cause an error if they exist
    for file in files_fail_ifexists:
        if str(file) in ignore_files:
            continue
        if pf(file).is_file():
            failed.append(f"File must be removed: {self._wrap_quotes(file)}")
        else:
            passed.append(f"File not found check: {self._wrap_quotes(file)}")

    # Files that cause a warning if they exist
    for file in files_warn_ifexists:
        if str(file) in ignore_files:
            continue
        if pf(file).is_file():
            warned.append(f"File should be removed: {self._wrap_quotes(file)}")
        else:
            passed.append(f"File not found check: {self._wrap_quotes(file)}")

    # Files that are ignored
    for file in ignore_files:
        ignored.append(f"File is ignored: {self._wrap_quotes(file)}")

    return {"passed": passed, "warned": warned, "failed": failed, "ignored": ignored}
