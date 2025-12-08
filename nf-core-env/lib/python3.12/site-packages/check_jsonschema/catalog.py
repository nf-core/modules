from __future__ import annotations

import typing as t


def _bitbucket_pipelines_url() -> str:
    return "https://api.bitbucket.org/schemas/pipelines-configuration"


def _githubusercontent_url(owner: str, repo: str, ref: str, path: str) -> str:
    return f"https://raw.githubusercontent.com/{owner}/{repo}/{ref}/{path}"


# this lists custom schemas which are *not* part of the catalog
CUSTOM_SCHEMA_NAMES = [
    "github-workflows-require-timeout",
]

# Known configs. The SchemaCatalog lists known schema URLs with their names.
# kept in alphabetical order by name
#
# Additional config could be associated with the schemas in the future.
SCHEMA_CATALOG: dict[str, dict[str, t.Any]] = {
    "azure-pipelines": {
        "url": _githubusercontent_url(
            "microsoft", "azure-pipelines-vscode", "main", "service-schema.json"
        ),
        "hook_config": {
            "name": "Validate Azure Pipelines",
            "description": (
                "Validate Azure Pipelines config against the schema provided "
                "by Microsoft"
            ),
            "add_args": [
                "--data-transform",
                "azure-pipelines",
                "--regex-variant",
                "nonunicode",
            ],
            "files": r"^(\.)?azure-pipelines\.(yml|yaml)$",
            "types": "yaml",
        },
    },
    "bamboo-spec": {
        "url": "https://json.schemastore.org/bamboo-spec.json",
        "hook_config": {
            "name": "Validate Bamboo Specs",
            "files": r"^bamboo-specs/.*\.(yml|yaml)$",
            "types": "yaml",
        },
    },
    "bitbucket-pipelines": {
        "url": _bitbucket_pipelines_url(),
        "hook_config": {
            "name": "Validate Bitbucket Pipelines",
            "files": r"bitbucket-pipelines\.(yml|yaml)$",
            "types": "yaml",
        },
    },
    "buildkite": {
        "url": _githubusercontent_url(
            "buildkite", "pipeline-schema", "main", "schema.json"
        ),
        "hook_config": {
            "name": "Validate Buildkite Pipelines",
            "description": (
                "Validate Buildkite Pipelines against the schema provided by Buildkite"
            ),
            "files": [
                r"buildkite\.(yml|yaml|json)",
                r"buildkite\.(.+)\.(yml|yaml|json)",
                r"(.*/)?\.buildkite/pipeline\.(yml|yaml|json)",
                r"(.*/)?\.buildkite/pipeline\.(.+)\.(yml|yaml|json)",
            ],
            "types_or": ["json", "yaml"],
        },
    },
    "circle-ci": {
        "url": _githubusercontent_url(
            "CircleCI-Public", "circleci-yaml-language-server", "main", "schema.json"
        ),
        "hook_config": {
            "name": "Validate CircleCI config",
            "description": (
                "Validate CircleCI config against the schema provided by SchemaStore"
            ),
            "files": r"^\.circleci/config\.(yml|yaml)$",
            "type": "yaml",
        },
    },
    "citation-file-format": {
        "url": _githubusercontent_url(
            "citation-file-format",
            "citation-file-format",
            "main",
            "schema.json",
        ),
        "hook_config": {
            "name": "Validate Citation File Format",
            "description": "Validate Citation File Format",
            "files": r"^CITATION.cff$",
        },
    },
    "cloudbuild": {
        "url": "https://json.schemastore.org/cloudbuild.json",
        "hook_config": {
            "name": "Validate Google Cloud Build config",
            "description": (
                "Validate Google Cloud Build config against the schema provided "
                "by SchemaStore"
            ),
            "files": r"^cloudbuild\.(yml|yaml|json)$",
            "types_or": ["json", "yaml"],
        },
    },
    "codecov": {
        "url": "https://www.schemastore.org/codecov.json",
        "hook_config": {
            "name": "Validate Codecov config",
            "files": [
                r"^((\.github|dev)/)?\.?codecov\.ya?ml$",
            ],
            "types": "yaml",
        },
    },
    "compose-spec": {
        "url": _githubusercontent_url(
            "compose-spec",
            "compose-spec",
            "master",
            "schema/compose-spec.json",
        ),
        "hook_config": {
            "name": "Validate Docker Compose files",
            "files": [
                r"([^/]*/)*docker-compose(\.[\.a-zA-Z0-9_-]*)*\.(yml|yaml)",
                r"([^/]*/)*compose(\.[\.a-zA-Z0-9_-]*)*\.(yml|yaml)",
            ],
            "types": "yaml",
        },
    },
    "dependabot": {
        "url": "https://json.schemastore.org/dependabot-2.0.json",
        "hook_config": {
            "name": "Validate Dependabot Config (v2)",
            "files": r"^\.github/dependabot\.(yml|yaml)$",
            "types": "yaml",
        },
    },
    "drone-ci": {
        "url": "https://json.schemastore.org/drone.json",
        "hook_config": {
            "name": "Validate Drone-CI Config",
            "files": r"^\.drone\.yml$",
            "types": "yaml",
        },
    },
    "github-actions": {
        "url": "https://json.schemastore.org/github-action",
        "hook_config": {
            "name": "Validate GitHub Actions",
            "files": [
                r"action\.(yml|yaml)",
                r"\.github/actions/(.+/)?action\.(yml|yaml)",
            ],
            "types": "yaml",
        },
    },
    "github-issue-config": {
        "url": "https://www.schemastore.org/github-issue-config.json",
        "hook_config": {
            "name": "Validate GitHub issue config",
            "files": [
                r"^\.github/ISSUE_TEMPLATE/config\.yml$",
            ],
            "types": "yaml",
        },
    },
    "github-issue-forms": {
        "url": "https://www.schemastore.org/github-issue-forms.json",
        "hook_config": {
            "name": "Validate GitHub issue forms",
            "files": [
                r"^\.github/ISSUE_TEMPLATE/(?!config\.ya?ml$).+$",
            ],
            "types": "yaml",
        },
    },
    "github-workflows": {
        "url": "https://json.schemastore.org/github-workflow",
        "hook_config": {
            "name": "Validate GitHub Workflows",
            "files": r"^\.github/workflows/[^/]+$",
            "types": "yaml",
        },
    },
    "gitlab-ci": {
        "url": (
            "https://gitlab.com/gitlab-org/gitlab/-/raw/master/app/assets/javascripts"
            "/editor/schema/ci.json"
        ),
        "hook_config": {
            "name": "Validate GitLab CI config",
            "add_args": [
                "--data-transform",
                "gitlab-ci",
                "--regex-variant",
                "nonunicode",
            ],
            "files": r"^.*\.gitlab-ci\.(yml|yaml)$",
            "types": "yaml",
        },
    },
    "meltano": {
        "url": _githubusercontent_url(
            "meltano",
            "meltano",
            "main",
            "src/meltano/schemas/meltano.schema.json",
        ),
        "hook_config": {
            "name": "Validate Meltano Config",
            "description": (
                "Validate Meltano config against the schema provided by Meltano"
            ),
            "files": [
                r".*meltano\.yml",
                r"meltano-manifest\.json",
                r"meltano-manifest\..+\.json",
            ],
            "types_or": ["json", "yaml"],
        },
    },
    "mergify": {
        "url": "https://docs.mergify.com/mergify-configuration-schema.json",
        "hook_config": {
            "name": "Validate Mergify config",
            "files": [
                r"\.mergify\.yml",
                r"\.mergify/config\.yml",
                r"\.github/mergify\.yml",
            ],
            "types": "yaml",
        },
    },
    "readthedocs": {
        "url": _githubusercontent_url(
            "readthedocs",
            "readthedocs.org",
            "main",
            "readthedocs/rtd_tests/fixtures/spec/v2/schema.json",
        ),
        "hook_config": {
            "name": "Validate ReadTheDocs Config",
            "description": (
                "Validate ReadTheDocs config against the schema provided by ReadTheDocs"
            ),
            "files": r"^\.readthedocs\.(yml|yaml)$",
            "types": "yaml",
        },
    },
    "renovate": {
        "url": "https://docs.renovatebot.com/renovate-schema.json",
        "hook_config": {
            "name": "Validate Renovate Config",
            "description": (
                "Validate Renovate config against the schema provided by "
                "Renovate (does not support renovate config in package.json)"
            ),
            "add_args": ["--regex-variant", "nonunicode"],
            "files": [
                r"renovate\.(json|json5)",
                r"\.(github|gitlab)/renovate\.(json|json5)",
                r"\.renovaterc",
                r"\.renovaterc\.(json|json5)",
            ],
        },
    },
    "snapcraft": {
        "url": _githubusercontent_url(
            "canonical", "snapcraft", "main", "schema/snapcraft.json"
        ),
        "hook_config": {
            "name": "Validate snapcraft files",
            "description": (
                "Validate snapcraft files against the schema provided by Canonical"
            ),
            "files": r"^(.+/)?snapcraft\.yaml$",
            "types": "yaml",
        },
    },
    "taskfile": {
        "url": "https://taskfile.dev/schema.json",
        "hook_config": {
            "name": "Validate Taskfile Config",
            "description": (
                "Validate Taskfile config against the schema provided by Task"
            ),
            "files": [
                r"Taskfile\.(yml|yaml)",
                r"taskfile\.(yml|yaml)",
                r"Taskfile\.dist\.(yml|yaml)",
                r"taskfile\.dist\.(yml|yaml)",
            ],
            "types": "yaml",
        },
    },
    "travis": {
        "url": "https://json.schemastore.org/travis",
        "hook_config": {
            "name": "Validate Travis Config",
            "files": r"^\.travis\.(yml|yaml)$",
            "types": "yaml",
        },
    },
    "woodpecker-ci": {
        "url": (
            "https://raw.githubusercontent.com/woodpecker-ci/woodpecker/main/pipeline"
            "/frontend/yaml/linter/schema/schema.json"
        ),
        "hook_config": {
            "name": "Validate Woodpecker Config",
            "files": [
                r"^\.woodpecker\.(yml|yaml)$",
                r"^\.woodpecker/.+\.(yml|yaml)$",
            ],
            "types": "yaml",
        },
    },
}
