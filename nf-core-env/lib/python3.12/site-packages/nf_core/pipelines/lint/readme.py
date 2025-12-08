import re
from pathlib import Path


def readme(self):
    """Repository ``README.md`` tests

    The ``README.md`` files for a project are very important and must meet some requirements:

    * Nextflow badge

      * If no Nextflow badge is found, a warning is given
      * If a badge is found but the version doesn't match the minimum version in the config file, the test fails
      * Example badge code:

        .. code-block:: md

           [![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.6-brightgreen.svg)](https://www.nextflow.io/)

    .. note:: This badge are a markdown image ``![alt-text](<image URL>)`` *inside* a markdown link ``[markdown image](<link URL>)``, so a bit fiddly to write.

    * Zenodo release

        * If pipeline is released but still contains a 'zenodo.XXXXXXX' tag, the test fails

    To disable this test, add the following to the pipeline's ``.nf-core.yml`` file:

    .. code-block:: yaml

        lint:
            readme: False

    To disable subsets of these tests, add the following to the pipeline's ``.nf-core.yml`` file:

    .. code-block:: yaml

            lint:
                readme:
                    - nextflow_badge
                    - nfcore_template_badge
                    - zenodo_release

    """
    passed = []
    warned = []
    failed = []

    # Remove field that should be ignored according to the linting config
    ignore_configs = self.lint_config.get("readme", []) if self.lint_config is not None else []

    with open(Path(self.wf_path, "README.md")) as fh:
        content = fh.read()

    if "nextflow_badge" not in ignore_configs:
        # Check that there is a readme badge showing the minimum required version of Nextflow
        # [![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
        # and that it has the correct version
        nf_badge_re = r"\[!\[Nextflow\]\(https://img\.shields\.io/badge/version-!?(?:%E2%89%A5|%3E%3D)([\d\.]+)-green\?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow\.io\)\]\(https://www\.nextflow\.io/\)"
        match = re.search(nf_badge_re, content)
        if match:
            nf_badge_version = match.group(1).strip("'\"")
            try:
                if nf_badge_version != self.minNextflowVersion:
                    raise AssertionError()
            except (AssertionError, KeyError):
                failed.append(
                    f"README Nextflow minimum version badge does not match config. Badge: `{nf_badge_version}`, "
                    f"Config: `{self.minNextflowVersion}`"
                )
            else:
                passed.append(
                    f"README Nextflow minimum version badge matched config. Badge: `{nf_badge_version}`, "
                    f"Config: `{self.minNextflowVersion}`"
                )
        else:
            warned.append("README did not have a Nextflow minimum version badge.")

    if "nfcore_template_badge" not in ignore_configs:
        # Check that there is a readme badge showing the current nf-core/tools template version
        # [![nf-core template version](https://img.shields.io/badge/nf--core_template-3.2.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.2.0)
        t_badge_re = r"\[!\[nf-core template version\]\(https://img\.shields\.io/badge/nf--core_template-([\d\.]+(\.dev[\d\.])?)-green\?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co\.re\)\]\(https://github\.com/nf-core/tools/releases/tag/([\d\.]+(\.dev[\d\.])?)\)"
        match = re.search(t_badge_re, content)
        if match:
            passed.append("README nf-core template version badge found.")
        else:
            warned.append("README did not have an nf-core template version badge.")

    if "zenodo_doi" not in ignore_configs:
        # Check that zenodo.XXXXXXX has been replaced with the zendo.DOI
        zenodo_re = r"/zenodo\.X+"
        match = re.search(zenodo_re, content)
        if match:
            warned.append(
                "README contains the placeholder `zenodo.XXXXXXX`. "
                "This should be replaced with the zenodo doi (after the first release)."
            )
        else:
            passed.append("README Zenodo placeholder was replaced with DOI.")

    return {"passed": passed, "warned": warned, "failed": failed}
