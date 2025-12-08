import os
from pathlib import Path

from ruamel.yaml import YAML


def version_consistency(self):
    """Pipeline and container version number consistency.

    .. note:: This test only runs when the ``--release`` flag is set for ``nf-core pipelines lint``,
              or ``$GITHUB_REF`` is equal to ``main``.

    This lint fetches the pipeline version number from four possible locations:

    * The pipeline config, ``manifest.version``
    * The docker container in the pipeline config, ``process.container``

        Some pipelines may not have this set on a pipeline level. If it is not found, it is ignored.

    * ``$GITHUB_REF``, if it looks like a release tag (``refs/tags/<something>``)
    * The YAML file .nf-core.yml

    The test then checks that:

    * The container name has a tag specified (eg. ``nfcore/pipeline:version``)
    * The pipeline version number is numeric (contains only numbers and dots)
    * That the version numbers all match one another
    """
    passed = []
    failed = []

    # Get the version definitions
    # Get version from nextflow.config
    versions = {}
    versions["manifest.version"] = self.nf_config.get("manifest.version", "").strip(" '\"")

    # Get version from the docker tag
    if self.nf_config.get("process.container", "") and ":" not in self.nf_config.get("process.container", ""):
        failed.append(f"Docker slug seems not to have a version tag: {self.nf_config.get('process.container', '')}")

    # Get config container tag (if set; one container per workflow)
    if self.nf_config.get("process.container", ""):
        versions["process.container"] = self.nf_config.get("process.container", "").strip(" '\"").split(":")[-1]

    # Get version from the $GITHUB_REF env var if this is a release
    if (
        os.environ.get("GITHUB_REF", "").startswith("refs/tags/")
        and os.environ.get("GITHUB_REPOSITORY", "") != "nf-core/tools"
    ):
        versions["GITHUB_REF"] = os.path.basename(os.environ["GITHUB_REF"].strip(" '\""))

    # Get version from the .nf-core.yml template
    yaml = YAML()
    nfcore_yml = yaml.load(Path(self.wf_path) / ".nf-core.yml")
    if nfcore_yml["template"] and "version" in nfcore_yml["template"]:
        versions["nfcore_yml.version"] = nfcore_yml["template"]["version"].strip(" '\"")

    # Check if they are all numeric
    for v_type, version in versions.items():
        if not version.replace(".", "").isdigit():
            failed.append(f"{v_type} was not numeric: {version}!")

    # Check if they are consistent
    if len(set(versions.values())) != 1:
        failed.append(
            "The versioning is not consistent between container, release tag and config. Found {}".format(
                ", ".join([f"{k} = {v}" for k, v in versions.items()])
            )
        )
    else:
        passed.append("Version tags are consistent: {}".format(", ".join([f"{k} = {v}" for k, v in versions.items()])))

    return {"passed": passed, "failed": failed}
