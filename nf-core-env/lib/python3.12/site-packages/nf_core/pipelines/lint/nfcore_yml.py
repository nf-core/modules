from pathlib import Path

from ruamel.yaml import YAML

from nf_core import __version__

REPOSITORY_TYPES = ["pipeline", "modules"]


def nfcore_yml(self) -> dict[str, list[str]]:
    """Repository ``.nf-core.yml`` tests

    The ``.nf-core.yml`` contains metadata for nf-core tools to correctly apply its features.

    * repository type:

        Check that the repository type is set.

    * nf core version:

        Check if the nf-core version is set to the latest version.

    """
    passed: list[str] = []
    warned: list[str] = []
    failed: list[str] = []
    ignored: list[str] = []

    yaml = YAML()

    # Remove field that should be ignored according to the linting config
    ignore_configs = self.lint_config.get(".nf-core", []) if self.lint_config is not None else []
    for ext in (".yml", ".yaml"):
        try:
            nf_core_yml = yaml.load(Path(self.wf_path) / f".nf-core{ext}")
            break
        except FileNotFoundError:
            continue
    else:
        raise FileNotFoundError("No `.nf-core.yml` file found.")

    if "repository_type" not in ignore_configs:
        # Check that the repository type is set in the .nf-core.yml
        if "repository_type" in nf_core_yml:
            repo_type = nf_core_yml["repository_type"]
            if repo_type not in REPOSITORY_TYPES:
                failed.append(
                    f"Repository type in `.nf-core.yml` is not valid. "
                    f"Should be one of `[{', '.join(REPOSITORY_TYPES)}]` but was `{repo_type}`"
                )
            else:
                passed.append(f"Repository type in `.nf-core.yml` is valid: `{repo_type}`")
        else:
            warned.append("Repository type not set in `.nf-core.yml`")
    else:
        ignored.append("`.nf-core.yml` variable ignored 'repository_type'")

    if "nf_core_version" not in ignore_configs:
        # Check that the nf-core version is set in the .nf-core.yml
        if "nf_core_version" in nf_core_yml:
            nf_core_version = nf_core_yml["nf_core_version"]
            if nf_core_version != __version__ and "dev" not in nf_core_version:
                warned.append(
                    f"nf-core version in `.nf-core.yml` is not set to the latest version. "
                    f"Should be `{__version__}` but was `{nf_core_version}`"
                )
            else:
                passed.append(f"nf-core version in `.nf-core.yml` is set to the latest version: `{nf_core_version}`")
        else:
            warned.append("nf-core version not set in `.nf-core.yml`")
    else:
        ignored.append("`.nf-core.yml` variable ignored 'nf_core_version'")

    return {"passed": passed, "warned": warned, "failed": failed, "ignored": ignored}
