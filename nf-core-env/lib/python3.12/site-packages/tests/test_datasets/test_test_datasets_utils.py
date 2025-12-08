import os
import unittest

import requests
import responses

from nf_core.test_datasets.test_datasets_utils import (
    MODULES_BRANCH_NAME,
    GithubApiEndpoints,
    create_download_url,
    create_pretty_nf_path,
    get_remote_branch_names,
    get_remote_tree_for_branch,
    list_files_by_branch,
)


def mock_get_pipelines_request(rsps: responses.RequestsMock) -> None:
    """
    Mock request to list pipelines
    """
    url = "https://raw.githubusercontent.com/nf-core/website/refs/heads/main/public/pipeline_names.json"
    resp_json = {
        "pipeline": [
            "airrflow",
            "ampliseq",
            "atacseq",
            "bacass",
            "bactmap",
        ]
    }

    # add dummy json response at given url for get request
    rsps.add(method="GET", url=url, json=resp_json, status=200)


def mock_gh_api_request(rsps: responses.RequestsMock, branch="mag") -> None:
    """
    Mock request to list files in a github branch
    """
    gh_urls = GithubApiEndpoints(gh_orga="nf-core", gh_repo="test-datasets")

    # output from requets for mag files
    resp_json = {
        "sha": "bd061d2421282afa00ae1c83b43151fda0e046b7",
        "url": "https://api.github.com/repos/nf-core/test-datasets/git/trees/bd061d2421282afa00ae1c83b43151fda0e046b7",
        "tree": [
            {
                "path": "README.md",
                "type": "blob",
                "url": "https://api.github.com/repos/nf-core/test-datasets/git/blobs/1ebf7eb543cbb09b875a8a1ef50e107c5cfa8310",
            },
            {
                "path": "assemblies",
                "type": "tree",
                "url": "https://api.github.com/repos/nf-core/test-datasets/git/trees/fdf5797e56e51b5e0bbd583e8499ef660204a194",
            },
            {
                "path": "assemblies/MEGAHIT-test_minigut.contigs.fa.gz",
                "type": "blob",
                "url": "https://api.github.com/repos/nf-core/test-datasets/git/blobs/83185a06b64158a5ea144b8a5a4c26499cd3f585",
            },
            {
                "path": "samplesheets/assembly_samplesheet.csv",
                "type": "blob",
                "url": "https://api.github.com/repos/nf-core/test-datasets/git/blobs/0b86b23c3f901103a1a72ec9953a67623cd90e3e",
            },
        ],
    }

    # Add dummy json response at gh url with ok status code
    rsps.add(method="GET", url=gh_urls.get_remote_tree_url_for_branch(branch), json=resp_json, status=200)


def mock_gh_api_download(rsps: responses.RequestsMock, branch: str, file: str) -> None:
    """
    Mock download request
    """
    gh_urls = GithubApiEndpoints(gh_orga="nf-core", gh_repo="test-datasets")
    url = gh_urls.get_file_download_url(branch, file)
    content = b"id,group,assembler,fasta\ntest_minigut,0,MEGAHIT,https://github.com/nf-core/test-datasets/raw/mag/assemblies/MEGAHIT-test_minigut.contigs.fa.gz\ntest_minigut,0,SPAdes,https://github.com/nf-core/test-datasets/raw/mag/assemblies/SPAdes-test_minigut_contigs.fasta.gz\ntest_minigut_sample2,0,MEGAHIT,https://github.com/nf-core/test-datasets/raw/mag/assemblies/MEGAHIT-test_minigut_sample2.contigs.fa.gz\ntest_minigut_sample2,0,SPAdes,https://github.com/nf-core/test-datasets/raw/mag/assemblies/SPAdes-test_minigut_sample2_contigs.fasta.gz\n"
    rsps.add(url=url, method="GET", body=content, content_type="text/plain")


class TestTestDatasetsUtils(unittest.TestCase):
    """Class for components tests"""

    def setUp(self):
        self.gh_urls = GithubApiEndpoints(gh_orga="nf-core", gh_repo="test-datasets")

    def test_modules_branch_name_changed(self):
        self.assertEqual(MODULES_BRANCH_NAME, "modules")

    def test_modules_branch_exists(self):
        url = self.gh_urls.get_remote_tree_url_for_branch(MODULES_BRANCH_NAME)
        resp = self._request_with_token(url)
        self.assertTrue(resp.ok)
        self.assertTrue(len(resp.json()) != 0)

    def test_get_branch_names(self):
        with responses.RequestsMock() as rsps:
            mock_get_pipelines_request(rsps)
            branch_names = get_remote_branch_names()
            self.assertTrue(len(branch_names) != 0)

    def test_get_remote_tree_for_branch(self):
        with responses.RequestsMock() as rsps:
            branch = "modules"
            mock_gh_api_request(rsps, branch)
            file_list = get_remote_tree_for_branch(branch)
            self.assertTrue(len(file_list) != 0)

    def test_list_files_by_branch(self):
        with responses.RequestsMock() as rsps:
            branch = "modules"
            mock_get_pipelines_request(rsps)  # Mocks fetching branch names
            mock_gh_api_request(rsps, branch)  # Mocks fetching file tree
            tree = list_files_by_branch(branch)
            self.assertTrue(len(tree.values()) != 0)
            self.assertTrue(tree.get(branch, None) is not None)

    def test_create_pretty_nf_path(self):
        nf_line_modules = create_pretty_nf_path("/path/to/file.xyz", is_module_dataset=True)
        nf_line_pipelines = create_pretty_nf_path("/path/to/file.xyz", is_module_dataset=False)
        self.assertEqual('params.modules_testdata_base_path + "/path/to/file.xyz"', nf_line_modules)
        self.assertEqual('params.pipelines_testdata_base_path + "/path/to/file.xyz"', nf_line_pipelines)

    def test_file_download(self):
        with responses.RequestsMock() as rsps:
            # create_download_url
            test_data_branch = "mag"
            small_test_file = "samplesheets/assembly_samplesheet.csv"
            mock_gh_api_download(rsps, test_data_branch, small_test_file)
            url = create_download_url(test_data_branch, small_test_file)
            resp = requests.get(url)
            self.assertTrue(resp.ok)
            self.assertTrue(len(resp.text) != 0)

    def test_github_endpoints(self):
        url_1 = self.gh_urls.get_remote_tree_url_for_branch(MODULES_BRANCH_NAME)
        url_2 = self.gh_urls.get_pipelines_list_url()
        url_3 = self.gh_urls.get_file_download_url(MODULES_BRANCH_NAME, "foo/bar/baz/qerljerlkmf")
        self.assertTrue(url_1 is not None)
        self.assertTrue(url_2 is not None)
        self.assertTrue(url_3 is not None)

        resp_1 = self._request_with_token(url_1)
        resp_2 = requests.get(url_2)
        resp_3 = requests.get(url_3)
        self.assertTrue(resp_1.ok)
        self.assertTrue(resp_2.ok)
        self.assertFalse(resp_3.ok)

    def _request_with_token(self, url):
        """Make a request with a GitHub token if available to mitigate issues with API rate limits."""
        headers = {}
        github_token = os.environ.get("GITHUB_TOKEN")
        if github_token:
            headers["authorization"] = f"Bearer {github_token}"

        resp = requests.get(url, headers=headers)
        return resp
