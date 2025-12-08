"""
Helper functions for tests
"""

import filecmp
import functools
import tempfile
from collections.abc import Callable
from pathlib import Path
from typing import Any

import responses
import yaml

import nf_core.modules
import nf_core.pipelines.create.create
from nf_core import __version__
from nf_core.utils import NFCoreTemplateConfig, NFCoreYamlConfig, custom_yaml_dumper

TEST_DATA_DIR = Path(__file__).parent / "data"
OLD_TRIMGALORE_SHA = "9b7a3bdefeaad5d42324aa7dd50f87bea1b04386"
OLD_TRIMGALORE_BRANCH = "mimic-old-trimgalore"
GITLAB_URL = "https://gitlab.com/nf-core/modules-test.git"
CROSS_ORGANIZATION_URL = "https://github.com/nf-core-test/modules.git"
CROSS_ORGANIZATION_BRANCH = "main"
GITLAB_REPO = "nf-core-test"
GITLAB_DEFAULT_BRANCH = "main"
GITLAB_SUBWORKFLOWS_BRANCH = "subworkflows"
GITLAB_SUBWORKFLOWS_ORG_PATH_BRANCH = "subworkflows-org-path"
OLD_SUBWORKFLOWS_SHA = "f3c078809a2513f1c95de14f6633fe1f03572fdb"
# Branch test stuff
GITLAB_BRANCH_TEST_BRANCH = "branch-tester"
GITLAB_BRANCH_ORG_PATH_BRANCH = "org-path"
GITLAB_BRANCH_TEST_OLD_SHA = "e772abc22c1ff26afdf377845c323172fb3c19ca"
GITLAB_BRANCH_TEST_NEW_SHA = "7d73e21f30041297ea44367f2b4fd4e045c0b991"
GITLAB_NFTEST_BRANCH = "nf-test-tests"


def with_temporary_folder(func: Callable[..., Any]) -> Callable[..., Any]:
    """
    Call the decorated function under the tempfile.TemporaryDirectory
    context manager. Pass the temporary directory name to the decorated
    function
    """

    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        with tempfile.TemporaryDirectory() as tmpdirname:
            return func(*args, tmpdirname, **kwargs)

    return wrapper


def with_temporary_file(func: Callable[..., Any]) -> Callable[..., Any]:
    """
    Call the decorated function under the tempfile.NamedTemporaryFile
    context manager. Pass the opened file handle to the decorated function
    """

    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        with tempfile.NamedTemporaryFile() as tmpfile:
            return func(*args, tmpfile, **kwargs)

    return wrapper


def mock_anaconda_api_calls(rsps: responses.RequestsMock, module: str, version: str) -> None:
    """Mock anaconda api calls for module"""
    anaconda_api_url = f"https://api.anaconda.org/package/bioconda/{module}"
    anaconda_mock = {
        "latest_version": version.split("--")[0],
        "summary": "",
        "doc_url": "http://test",
        "dev_url": "http://test",
        "files": [{"version": version.split("--")[0]}],
        "license": "MIT",
    }
    rsps.get(anaconda_api_url, json=anaconda_mock, status=200)


def mock_biocontainers_api_calls(rsps: responses.RequestsMock, module: str, version: str) -> None:
    """Mock biocontainers api calls for module"""
    biocontainers_api_url = (
        f"https://api.biocontainers.pro/ga4gh/trs/v2/tools/{module}/versions/{module}-{version.split('--')[0]}"
    )
    biocontainers_mock = {
        "images": [
            {
                "image_type": "Singularity",
                "image_name": f"https://depot.galaxyproject.org/singularity/{module}:{version}",
                "updated": "2021-09-04T00:00:00Z",
            },
            {
                "image_type": "Docker",
                "image_name": f"biocontainers/{module}:{version}",
                "updated": "2021-09-04T00:00:00Z",
            },
        ],
    }
    rsps.get(biocontainers_api_url, json=biocontainers_mock, status=200)


def mock_biotools_api_calls(rsps: responses.RequestsMock, module: str) -> None:
    """Mock biotools api calls for module"""
    biotools_api_url = f"https://bio.tools/api/t/?q={module}&format=json"
    biotools_mock = {
        "list": [
            {
                "name": "Bpipe",
                "biotoolsCURIE": "biotools:bpipe",
                "function": [
                    {
                        "input": [
                            {
                                "data": {"uri": "http://edamontology.org/data_0848", "term": "Raw sequence"},
                                "format": [
                                    {"uri": "http://edamontology.org/format_2182", "term": "FASTQ-like format (text)"},
                                    {"uri": "http://edamontology.org/format_2573", "term": "SAM"},
                                ],
                            }
                        ],
                        "output": [
                            {
                                "data": {"uri": "http://edamontology.org/data_2955", "term": "Sequence report"},
                                "format": [{"uri": "http://edamontology.org/format_2331", "term": "HTML"}],
                            }
                        ],
                    }
                ],
            }
        ],
    }
    rsps.get(biotools_api_url, json=biotools_mock, status=200)


def create_tmp_pipeline(no_git: bool = False) -> tuple[Path, Path, str, Path]:
    """Create a new Pipeline for testing"""

    tmp_dir = Path(tempfile.TemporaryDirectory().name)
    root_repo_dir = Path(__file__).resolve().parent.parent
    template_dir = root_repo_dir / "nf_core" / "pipeline-template"
    pipeline_name = "testpipeline"
    pipeline_dir = tmp_dir / pipeline_name
    pipeline_dir.mkdir(parents=True)

    nf_core_yml = NFCoreYamlConfig(
        nf_core_version=__version__,
        repository_type="modules",
        org_path="nf-core",
        lint=None,
        template=NFCoreTemplateConfig(
            name="testpipeline",
            author="me",
            description="it is mine",
            org="nf-core",
            version=None,
            force=True,
            is_nfcore=None,
            skip_features=None,
            outdir=None,
        ),
        bump_version=None,
    )
    with open(str(Path(pipeline_dir, ".nf-core.yml")), "w") as fh:
        yaml.dump(nf_core_yml.model_dump(), fh, Dumper=custom_yaml_dumper())

    nf_core.pipelines.create.create.PipelineCreate(
        pipeline_name, "it is mine", "me", no_git=no_git, outdir=pipeline_dir, force=True
    ).init_pipeline()

    # return values to instance variables for later use in test methods
    return tmp_dir, template_dir, pipeline_name, pipeline_dir


def cmp_component(dir1: Path, dir2: Path) -> bool:
    """Compare two versions of the same component"""
    files = ["main.nf", "meta.yml"]
    return all(filecmp.cmp(dir1 / f, dir2 / f, shallow=False) for f in files)
