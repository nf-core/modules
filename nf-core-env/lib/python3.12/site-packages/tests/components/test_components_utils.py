import importlib
import os
from unittest import mock

import responses

import nf_core.components.components_utils

from ..test_components import TestComponents
from ..utils import mock_biotools_api_calls


class TestTestComponentsUtils(TestComponents):
    def test_get_biotools_id(self):
        """Test getting the bio.tools ID for a tool"""
        with responses.RequestsMock() as rsps:
            mock_biotools_api_calls(rsps, "bpipe")
            response = nf_core.components.components_utils.get_biotools_response("bpipe")
            id = nf_core.components.components_utils.get_biotools_id(response, "bpipe")
            assert id == "biotools:bpipe"

    def test_get_biotools_id_warn(self):
        """Test getting the bio.tools ID for a tool and failing"""
        with responses.RequestsMock() as rsps:
            mock_biotools_api_calls(rsps, "bpipe")
            response = nf_core.components.components_utils.get_biotools_response("bpipe")
            nf_core.components.components_utils.get_biotools_id(response, "test")
            assert "Could not find a bio.tools ID for 'test'" in self.caplog.text

    def test_get_biotools_ch_info(self):
        """Test getting the bio.tools channel information for a tool"""
        with responses.RequestsMock() as rsps:
            mock_biotools_api_calls(rsps, "bpipe")
            response = nf_core.components.components_utils.get_biotools_response("bpipe")
            inputs, outputs = nf_core.components.components_utils.get_channel_info_from_biotools(response, "bpipe")
            assert inputs == {
                "raw_sequence": (
                    [
                        "http://edamontology.org/data_0848",
                        "http://edamontology.org/format_2182",
                        "http://edamontology.org/format_2573",
                    ],
                    ["Raw sequence", "FASTQ-like format (text)", "SAM"],
                    ["fastq-like", "sam"],
                )
            }
            assert outputs == {
                "sequence_report": (
                    ["http://edamontology.org/data_2955", "http://edamontology.org/format_2331"],
                    ["Sequence report", "HTML"],
                    ["html"],
                )
            }

    def test_get_biotools_ch_info_warn(self):
        """Test getting the bio.tools channel information for a tool and failing"""
        with responses.RequestsMock() as rsps:
            mock_biotools_api_calls(rsps, "bpipe")
            response = nf_core.components.components_utils.get_biotools_response("bpipe")
            nf_core.components.components_utils.get_channel_info_from_biotools(response, "test")
            assert "Could not find an EDAM ontology term for 'test'" in self.caplog.text

    def test_environment_variables_override(self):
        """Test environment variables override default values"""
        mock_env = {
            "NF_CORE_MODULES_NAME": "custom-name",
            "NF_CORE_MODULES_REMOTE": "https://custom-repo.git",
            "NF_CORE_MODULES_DEFAULT_BRANCH": "custom-branch",
        }

        try:
            with mock.patch.dict(os.environ, mock_env):
                importlib.reload(nf_core.components.constants)
                assert nf_core.components.constants.NF_CORE_MODULES_NAME == mock_env["NF_CORE_MODULES_NAME"]
                assert nf_core.components.constants.NF_CORE_MODULES_REMOTE == mock_env["NF_CORE_MODULES_REMOTE"]
                assert (
                    nf_core.components.constants.NF_CORE_MODULES_DEFAULT_BRANCH
                    == mock_env["NF_CORE_MODULES_DEFAULT_BRANCH"]
                )
        finally:
            importlib.reload(nf_core.components.constants)
