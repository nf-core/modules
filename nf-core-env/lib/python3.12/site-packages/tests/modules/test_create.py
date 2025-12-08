import os
import shutil
from pathlib import Path
from unittest import mock

import pytest
import requests_cache
import responses
import yaml
from git.repo import Repo

import nf_core.modules.create
from tests.utils import (
    GITLAB_SUBWORKFLOWS_ORG_PATH_BRANCH,
    GITLAB_URL,
    mock_anaconda_api_calls,
    mock_biocontainers_api_calls,
    mock_biotools_api_calls,
)

from ..test_modules import TestModules


class TestModulesCreate(TestModules):
    def test_modules_create_succeed(self):
        """Succeed at creating the TrimGalore! module"""
        with responses.RequestsMock() as rsps:
            mock_anaconda_api_calls(rsps, "trim-galore", "0.6.7")
            mock_biocontainers_api_calls(rsps, "trim-galore", "0.6.7")
            module_create = nf_core.modules.create.ModuleCreate(
                self.pipeline_dir, "trimgalore", "@author", "process_single", True, True, conda_name="trim-galore"
            )
            with requests_cache.disabled():
                module_create.create()
        assert os.path.exists(os.path.join(self.pipeline_dir, "modules", "local", "trimgalore/main.nf"))

    def test_modules_create_fail_exists(self):
        """Fail at creating the same module twice"""
        with responses.RequestsMock() as rsps:
            mock_anaconda_api_calls(rsps, "trim-galore", "0.6.7")
            mock_biocontainers_api_calls(rsps, "trim-galore", "0.6.7")
            module_create = nf_core.modules.create.ModuleCreate(
                self.pipeline_dir, "trimgalore", "@author", "process_single", False, False, conda_name="trim-galore"
            )
            with requests_cache.disabled():
                module_create.create()
            with pytest.raises(UserWarning) as excinfo:
                with requests_cache.disabled():
                    module_create.create()
        assert "module directory exists:" in str(excinfo.value)

    def test_modules_create_nfcore_modules(self):
        """Create a module in nf-core/modules clone"""
        with responses.RequestsMock() as rsps:
            mock_anaconda_api_calls(rsps, "fastqc", "0.11.9")
            mock_biocontainers_api_calls(rsps, "fastqc", "0.11.9")
            module_create = nf_core.modules.create.ModuleCreate(
                self.nfcore_modules, "fastqc", "@author", "process_low", False, False
            )
            with requests_cache.disabled():
                module_create.create()
        assert os.path.exists(os.path.join(self.nfcore_modules, "modules", "nf-core", "fastqc", "main.nf"))
        assert os.path.exists(
            os.path.join(self.nfcore_modules, "modules", "nf-core", "fastqc", "tests", "main.nf.test")
        )

    def test_modules_create_nfcore_modules_subtool(self):
        """Create a tool/subtool module in a nf-core/modules clone"""
        with responses.RequestsMock() as rsps:
            mock_anaconda_api_calls(rsps, "star", "2.8.10a")
            mock_biocontainers_api_calls(rsps, "star", "2.8.10a")
            module_create = nf_core.modules.create.ModuleCreate(
                self.nfcore_modules, "star/index", "@author", "process_medium", False, False
            )
            with requests_cache.disabled():
                module_create.create()
        assert os.path.exists(os.path.join(self.nfcore_modules, "modules", "nf-core", "star", "index", "main.nf"))
        assert os.path.exists(
            os.path.join(self.nfcore_modules, "modules", "nf-core", "star", "index", "tests", "main.nf.test")
        )

    @mock.patch("rich.prompt.Confirm.ask")
    def test_modules_migrate(self, mock_rich_ask):
        """Create a module with the --migrate-pytest option to convert pytest to nf-test"""
        pytest_dir = Path(self.nfcore_modules, "tests", "modules", "nf-core", "samtools", "sort")
        module_dir = Path(self.nfcore_modules, "modules", "nf-core", "samtools", "sort")

        # Clone modules repo with pytests
        shutil.rmtree(self.nfcore_modules)
        Repo.clone_from(GITLAB_URL, self.nfcore_modules, branch=GITLAB_SUBWORKFLOWS_ORG_PATH_BRANCH)
        with open(module_dir / "main.nf") as fh:
            old_main_nf = fh.read()
        with open(module_dir / "meta.yml") as fh:
            old_meta_yml = fh.read()

        # Create a module with --migrate-pytest
        mock_rich_ask.return_value = True
        module_create = nf_core.modules.create.ModuleCreate(self.nfcore_modules, "samtools/sort", migrate_pytest=True)
        module_create.create()

        with open(module_dir / "main.nf") as fh:
            new_main_nf = fh.read()
        with open(module_dir / "meta.yml") as fh:
            new_meta_yml = fh.read()
        nextflow_config = module_dir / "tests" / "nextflow.config"

        # Check that old files have been copied to the new module
        assert old_main_nf == new_main_nf
        assert old_meta_yml == new_meta_yml
        assert nextflow_config.is_file()

        # Check that pytest folder is deleted
        assert not pytest_dir.is_dir()

        # Check that pytest_modules.yml is updated
        with open(Path(self.nfcore_modules, "tests", "config", "pytest_modules.yml")) as fh:
            modules_yml = yaml.safe_load(fh)
        assert "samtools/sort" not in modules_yml.keys()

    @mock.patch("rich.prompt.Confirm.ask")
    def test_modules_migrate_no_delete(self, mock_rich_ask):
        """Create a module with the --migrate-pytest option to convert pytest to nf-test.
        Test that pytest directory is not deleted."""
        pytest_dir = Path(self.nfcore_modules, "tests", "modules", "nf-core", "samtools", "sort")

        # Clone modules repo with pytests
        shutil.rmtree(self.nfcore_modules)
        Repo.clone_from(GITLAB_URL, self.nfcore_modules, branch=GITLAB_SUBWORKFLOWS_ORG_PATH_BRANCH)

        # Create a module with --migrate-pytest
        mock_rich_ask.return_value = False
        module_create = nf_core.modules.create.ModuleCreate(self.nfcore_modules, "samtools/sort", migrate_pytest=True)
        module_create.create()

        # Check that pytest folder is not deleted
        assert pytest_dir.is_dir()

        # Check that pytest_modules.yml is updated
        with open(Path(self.nfcore_modules, "tests", "config", "pytest_modules.yml")) as fh:
            modules_yml = yaml.safe_load(fh)
        assert "samtools/sort" not in modules_yml.keys()

    @mock.patch("rich.prompt.Confirm.ask")
    def test_modules_migrate_symlink(self, mock_rich_ask):
        """Create a module with the --migrate-pytest option to convert pytest with symlinks to nf-test.
        Test that the symlink is deleted and the file is copied."""

        pytest_dir = Path(self.nfcore_modules, "tests", "modules", "nf-core", "samtools", "sort")
        module_dir = Path(self.nfcore_modules, "modules", "nf-core", "samtools", "sort")

        # Clone modules repo with pytests
        shutil.rmtree(self.nfcore_modules)
        Repo.clone_from(GITLAB_URL, self.nfcore_modules, branch=GITLAB_SUBWORKFLOWS_ORG_PATH_BRANCH)

        # Create a symlinked file in the pytest directory
        symlink_file = pytest_dir / "symlink_file.txt"
        symlink_file.symlink_to(module_dir / "main.nf")

        # Create a module with --migrate-pytest
        mock_rich_ask.return_value = True
        module_create = nf_core.modules.create.ModuleCreate(self.nfcore_modules, "samtools/sort", migrate_pytest=True)
        module_create.create()

        # Check that symlink is deleted
        assert not symlink_file.is_symlink()

    def test_modules_meta_yml_structure_biotools_meta(self):
        """Test the structure of the module meta.yml file when it was generated with INFORMATION from bio.tools and WITH a meta."""
        with responses.RequestsMock() as rsps:
            mock_anaconda_api_calls(rsps, "bpipe", "0.9.13--hdfd78af_0")
            mock_biocontainers_api_calls(rsps, "bpipe", "0.9.13--hdfd78af_0")
            mock_biotools_api_calls(rsps, "bpipe")
            module_create = nf_core.modules.create.ModuleCreate(
                self.nfcore_modules, "bpipe/test", "@author", "process_single", has_meta=True, force=True
            )
            module_create.create()

        expected_yml = {
            "name": "bpipe_test",
            "description": "write your description here",
            "keywords": ["sort", "example", "genomics"],
            "tools": [
                {
                    "bpipe": {
                        "description": "",
                        "homepage": "http://test",
                        "documentation": "http://test",
                        "tool_dev_url": "http://test",
                        "doi": "",
                        "licence": ["MIT"],
                        "identifier": "biotools:bpipe",
                    }
                }
            ],
            "input": [
                [
                    {
                        "meta": {
                            "type": "map",
                            "description": "Groovy Map containing sample information. e.g. `[ id:'sample1' ]`",
                        }
                    },
                    {
                        "raw_sequence": {
                            "type": "file",
                            "description": "raw_sequence file",
                            "pattern": "*.{fastq-like,sam}",
                            "ontologies": [
                                {"edam": "http://edamontology.org/data_0848"},
                                {"edam": "http://edamontology.org/format_2182"},
                                {"edam": "http://edamontology.org/format_2573"},
                            ],
                        }
                    },
                ]
            ],
            "output": {
                "sequence_report": [
                    [
                        {
                            "meta": {
                                "type": "map",
                                "description": "Groovy Map containing sample information. e.g. `[ id:'sample1' ]`",
                            }
                        },
                        {
                            "*.{html}": {
                                "type": "file",
                                "description": "sequence_report file",
                                "pattern": "*.{html}",
                                "ontologies": [
                                    {"edam": "http://edamontology.org/data_2955"},
                                    {"edam": "http://edamontology.org/format_2331"},
                                ],
                            }
                        },
                    ]
                ],
                "versions_bpipe": [
                    [
                        {"${task.process}": {"type": "string", "description": "The name of the process"}},
                        {"bpipe": {"type": "string", "description": "The name of the tool"}},
                        {
                            "bpipe --version": {
                                "type": "eval",
                                "description": "The expression to obtain the version of the tool",
                            }
                        },
                    ]
                ],
            },
            "topics": {
                "versions": [
                    [
                        {"${task.process}": {"type": "string", "description": "The name of the process"}},
                        {"bpipe": {"type": "string", "description": "The name of the tool"}},
                        {
                            "bpipe --version": {
                                "type": "eval",
                                "description": "The expression to obtain the version of the tool",
                            }
                        },
                    ]
                ],
            },
            "authors": ["@author"],
            "maintainers": ["@author"],
        }

        with open(Path(self.nfcore_modules, "modules", "nf-core", "bpipe", "test", "meta.yml")) as fh:
            meta_yml = yaml.safe_load(fh)

        assert meta_yml == expected_yml

    def test_modules_meta_yml_structure_biotools_nometa(self):
        """Test the structure of the module meta.yml file when it was generated with INFORMATION from bio.tools and WITHOUT a meta."""
        with responses.RequestsMock() as rsps:
            mock_anaconda_api_calls(rsps, "bpipe", "0.9.13--hdfd78af_0")
            mock_biocontainers_api_calls(rsps, "bpipe", "0.9.13--hdfd78af_0")
            mock_biotools_api_calls(rsps, "bpipe")
            module_create = nf_core.modules.create.ModuleCreate(
                self.nfcore_modules, "bpipe/test", "@author", "process_single", has_meta=False, force=True
            )
            module_create.create()

        expected_yml = {
            "name": "bpipe_test",
            "description": "write your description here",
            "keywords": ["sort", "example", "genomics"],
            "tools": [
                {
                    "bpipe": {
                        "description": "",
                        "homepage": "http://test",
                        "documentation": "http://test",
                        "tool_dev_url": "http://test",
                        "doi": "",
                        "licence": ["MIT"],
                        "identifier": "biotools:bpipe",
                    }
                }
            ],
            "input": [
                {
                    "raw_sequence": {
                        "type": "file",
                        "description": "raw_sequence file",
                        "pattern": "*.{fastq-like,sam}",
                        "ontologies": [
                            {"edam": "http://edamontology.org/data_0848"},
                            {"edam": "http://edamontology.org/format_2182"},
                            {"edam": "http://edamontology.org/format_2573"},
                        ],
                    }
                }
            ],
            "output": {
                "sequence_report": [
                    {
                        "*.{html}": {
                            "type": "file",
                            "description": "sequence_report file",
                            "pattern": "*.{html}",
                            "ontologies": [
                                {"edam": "http://edamontology.org/data_2955"},
                                {"edam": "http://edamontology.org/format_2331"},
                            ],
                        }
                    }
                ],
                "versions_bpipe": [
                    [
                        {"${task.process}": {"type": "string", "description": "The name of the process"}},
                        {"bpipe": {"type": "string", "description": "The name of the tool"}},
                        {
                            "bpipe --version": {
                                "type": "eval",
                                "description": "The expression to obtain the version of the tool",
                            }
                        },
                    ]
                ],
            },
            "topics": {
                "versions": [
                    [
                        {"${task.process}": {"type": "string", "description": "The name of the process"}},
                        {"bpipe": {"type": "string", "description": "The name of the tool"}},
                        {
                            "bpipe --version": {
                                "type": "eval",
                                "description": "The expression to obtain the version of the tool",
                            }
                        },
                    ]
                ],
            },
            "authors": ["@author"],
            "maintainers": ["@author"],
        }

        with open(Path(self.nfcore_modules, "modules", "nf-core", "bpipe", "test", "meta.yml")) as fh:
            meta_yml = yaml.safe_load(fh)

        assert meta_yml == expected_yml

    @mock.patch("nf_core.utils.anaconda_package")
    @mock.patch("nf_core.utils.get_biocontainer_tag")
    @mock.patch("nf_core.components.components_utils.get_biotools_response")
    @mock.patch("rich.prompt.Confirm.ask")
    def test_modules_meta_yml_structure_template_meta(
        self, mock_rich_ask, mock_biotools_response, mock_biocontainer_tag, mock_anaconda_package
    ):
        """Test the structure of the module meta.yml file when it was generated with TEMPLATE data and WITH a meta."""
        mock_biotools_response.return_value = {}
        mock_biocontainer_tag.return_value = {}
        mock_anaconda_package.return_value = {}
        mock_rich_ask.return_value = False  # Don't provide Bioconda package name
        module_create = nf_core.modules.create.ModuleCreate(
            self.nfcore_modules, "test", "@author", "process_single", has_meta=True, empty_template=False
        )
        module_create.create()

        expected_yml = {
            "name": "test",
            "description": "write your description here",
            "keywords": ["sort", "example", "genomics"],
            "tools": [
                {
                    "test": {
                        "description": "",
                        "homepage": "",
                        "documentation": "",
                        "tool_dev_url": "",
                        "doi": "",
                        "licence": None,
                        "identifier": None,
                    }
                }
            ],
            "input": [
                [
                    {
                        "meta": {
                            "type": "map",
                            "description": "Groovy Map containing sample information\ne.g. `[ id:'sample1' ]`\n",
                        }
                    },
                    {
                        "bam": {
                            "type": "file",
                            "description": "Sorted BAM/CRAM/SAM file",
                            "pattern": "*.{bam,cram,sam}",
                            "ontologies": [
                                {"edam": "http://edamontology.org/format_2572"},
                                {"edam": "http://edamontology.org/format_2573"},
                                {"edam": "http://edamontology.org/format_3462"},
                            ],
                        }
                    },
                ]
            ],
            "output": {
                "bam": [
                    [
                        {
                            "meta": {
                                "type": "map",
                                "description": "Groovy Map containing sample information\ne.g. `[ id:'sample1' ]`\n",
                            }
                        },
                        {
                            "*.bam": {
                                "type": "file",
                                "description": "Sorted BAM/CRAM/SAM file",
                                "pattern": "*.{bam,cram,sam}",
                                "ontologies": [
                                    {"edam": "http://edamontology.org/format_2572"},
                                    {"edam": "http://edamontology.org/format_2573"},
                                    {"edam": "http://edamontology.org/format_3462"},
                                ],
                            }
                        },
                    ]
                ],
                "versions_test": [
                    [
                        {"${task.process}": {"type": "string", "description": "The name of the process"}},
                        {"test": {"type": "string", "description": "The name of the tool"}},
                        {
                            "test --version": {
                                "type": "eval",
                                "description": "The expression to obtain the version of the tool",
                            }
                        },
                    ]
                ],
            },
            "topics": {
                "versions": [
                    [
                        {"${task.process}": {"type": "string", "description": "The name of the process"}},
                        {"test": {"type": "string", "description": "The name of the tool"}},
                        {
                            "test --version": {
                                "type": "eval",
                                "description": "The expression to obtain the version of the tool",
                            }
                        },
                    ]
                ],
            },
            "authors": ["@author"],
            "maintainers": ["@author"],
        }

        with open(Path(self.nfcore_modules, "modules", "nf-core", "test", "meta.yml")) as fh:
            meta_yml = yaml.safe_load(fh)

        assert meta_yml == expected_yml

    @mock.patch("nf_core.utils.anaconda_package")
    @mock.patch("nf_core.utils.get_biocontainer_tag")
    @mock.patch("nf_core.components.components_utils.get_biotools_response")
    @mock.patch("rich.prompt.Confirm.ask")
    def test_modules_meta_yml_structure_template_nometa(
        self, mock_rich_ask, mock_biotools_response, mock_biocontainer_tag, mock_anaconda_package
    ):
        """Test the structure of the module meta.yml file when it was generated with TEMPLATE data and WITHOUT a meta."""
        mock_biotools_response.return_value = {}
        mock_biocontainer_tag.return_value = {}
        mock_anaconda_package.return_value = {}
        mock_rich_ask.return_value = False  # Don't provide Bioconda package name
        module_create = nf_core.modules.create.ModuleCreate(
            self.nfcore_modules, "test", "@author", "process_single", has_meta=False, empty_template=False
        )
        module_create.create()

        expected_yml = {
            "name": "test",
            "description": "write your description here",
            "keywords": ["sort", "example", "genomics"],
            "tools": [
                {
                    "test": {
                        "description": "",
                        "homepage": "",
                        "documentation": "",
                        "tool_dev_url": "",
                        "doi": "",
                        "licence": None,
                        "identifier": None,
                    }
                }
            ],
            "input": [
                {
                    "bam": {
                        "type": "file",
                        "description": "Sorted BAM/CRAM/SAM file",
                        "pattern": "*.{bam,cram,sam}",
                        "ontologies": [
                            {"edam": "http://edamontology.org/format_2572"},
                            {"edam": "http://edamontology.org/format_2573"},
                            {"edam": "http://edamontology.org/format_3462"},
                        ],
                    }
                }
            ],
            "output": {
                "bam": [
                    {
                        "*.bam": {
                            "type": "file",
                            "description": "Sorted BAM/CRAM/SAM file",
                            "pattern": "*.{bam,cram,sam}",
                            "ontologies": [
                                {"edam": "http://edamontology.org/format_2572"},
                                {"edam": "http://edamontology.org/format_2573"},
                                {"edam": "http://edamontology.org/format_3462"},
                            ],
                        }
                    }
                ],
                "versions_test": [
                    [
                        {"${task.process}": {"type": "string", "description": "The name of the process"}},
                        {"test": {"type": "string", "description": "The name of the tool"}},
                        {
                            "test --version": {
                                "type": "eval",
                                "description": "The expression to obtain the version of the tool",
                            }
                        },
                    ]
                ],
            },
            "topics": {
                "versions": [
                    [
                        {"${task.process}": {"type": "string", "description": "The name of the process"}},
                        {"test": {"type": "string", "description": "The name of the tool"}},
                        {
                            "test --version": {
                                "type": "eval",
                                "description": "The expression to obtain the version of the tool",
                            }
                        },
                    ]
                ],
            },
            "authors": ["@author"],
            "maintainers": ["@author"],
        }

        with open(Path(self.nfcore_modules, "modules", "nf-core", "test", "meta.yml")) as fh:
            meta_yml = yaml.safe_load(fh)

        assert meta_yml == expected_yml

    @mock.patch("nf_core.utils.anaconda_package")
    @mock.patch("nf_core.utils.get_biocontainer_tag")
    @mock.patch("nf_core.components.components_utils.get_biotools_response")
    @mock.patch("rich.prompt.Confirm.ask")
    def test_modules_meta_yml_structure_empty_meta(
        self, mock_rich_ask, mock_biotools_response, mock_biocontainer_tag, mock_anaconda_package
    ):
        """Test the structure of the module meta.yml file when it was generated with an EMPTY template and WITH a meta."""
        mock_biotools_response.return_value = {}
        mock_biocontainer_tag.return_value = {}
        mock_anaconda_package.return_value = {}
        mock_rich_ask.return_value = False  # Don't provide Bioconda package name
        module_create = nf_core.modules.create.ModuleCreate(
            self.nfcore_modules, "test", "@author", "process_single", has_meta=True, empty_template=True
        )
        module_create.create()

        expected_yml = {
            "name": "test",
            "description": "write your description here",
            "keywords": ["sort", "example", "genomics"],
            "tools": [
                {
                    "test": {
                        "description": "",
                        "homepage": "",
                        "documentation": "",
                        "tool_dev_url": "",
                        "doi": "",
                        "licence": None,
                        "identifier": None,
                    }
                }
            ],
            "input": [
                [
                    {
                        "meta": {
                            "type": "map",
                            "description": "Groovy Map containing sample information. e.g. `[ id:'sample1' ]`",
                        }
                    },
                    {"input": {"type": "file", "description": "", "pattern": "", "ontologies": [{"edam": ""}]}},
                ]
            ],
            "output": {
                "output": [
                    [
                        {
                            "meta": {
                                "type": "map",
                                "description": "Groovy Map containing sample information. e.g. `[ id:'sample1' ]`",
                            }
                        },
                        {"*": {"type": "file", "description": "", "pattern": "", "ontologies": [{"edam": ""}]}},
                    ]
                ],
                "versions_test": [
                    [
                        {"${task.process}": {"type": "string", "description": "The name of the process"}},
                        {"test": {"type": "string", "description": "The name of the tool"}},
                        {
                            "test --version": {
                                "type": "eval",
                                "description": "The expression to obtain the version of the tool",
                            }
                        },
                    ]
                ],
            },
            "topics": {
                "versions": [
                    [
                        {"${task.process}": {"type": "string", "description": "The name of the process"}},
                        {"test": {"type": "string", "description": "The name of the tool"}},
                        {
                            "test --version": {
                                "type": "eval",
                                "description": "The expression to obtain the version of the tool",
                            }
                        },
                    ]
                ],
            },
            "authors": ["@author"],
            "maintainers": ["@author"],
        }

        with open(Path(self.nfcore_modules, "modules", "nf-core", "test", "meta.yml")) as fh:
            meta_yml = yaml.safe_load(fh)

        assert meta_yml == expected_yml

    @mock.patch("nf_core.utils.anaconda_package")
    @mock.patch("nf_core.utils.get_biocontainer_tag")
    @mock.patch("nf_core.components.components_utils.get_biotools_response")
    @mock.patch("rich.prompt.Confirm.ask")
    def test_modules_meta_yml_structure_empty_nometa(
        self, mock_rich_ask, mock_biotools_response, mock_biocontainer_tag, mock_anaconda_package
    ):
        """Test the structure of the module meta.yml file when it was generated with an EMPTY template and WITHOUT a meta."""
        mock_biotools_response.return_value = {}
        mock_biocontainer_tag.return_value = {}
        mock_anaconda_package.return_value = {}
        mock_rich_ask.return_value = False  # Don't provide Bioconda package name
        module_create = nf_core.modules.create.ModuleCreate(
            self.nfcore_modules, "test", "@author", "process_single", has_meta=False, empty_template=True
        )
        module_create.create()

        expected_yml = {
            "name": "test",
            "description": "write your description here",
            "keywords": ["sort", "example", "genomics"],
            "tools": [
                {
                    "test": {
                        "description": "",
                        "homepage": "",
                        "documentation": "",
                        "tool_dev_url": "",
                        "doi": "",
                        "licence": None,
                        "identifier": None,
                    }
                }
            ],
            "input": [{"input": {"type": "file", "description": "", "pattern": "", "ontologies": [{"edam": ""}]}}],
            "output": {
                "output": [{"*": {"type": "file", "description": "", "pattern": "", "ontologies": [{"edam": ""}]}}],
                "versions_test": [
                    [
                        {"${task.process}": {"type": "string", "description": "The name of the process"}},
                        {"test": {"type": "string", "description": "The name of the tool"}},
                        {
                            "test --version": {
                                "type": "eval",
                                "description": "The expression to obtain the version of the tool",
                            }
                        },
                    ]
                ],
            },
            "topics": {
                "versions": [
                    [
                        {"${task.process}": {"type": "string", "description": "The name of the process"}},
                        {"test": {"type": "string", "description": "The name of the tool"}},
                        {
                            "test --version": {
                                "type": "eval",
                                "description": "The expression to obtain the version of the tool",
                            }
                        },
                    ]
                ],
            },
            "authors": ["@author"],
            "maintainers": ["@author"],
        }

        with open(Path(self.nfcore_modules, "modules", "nf-core", "test", "meta.yml")) as fh:
            meta_yml = yaml.safe_load(fh)

        assert meta_yml == expected_yml
