from pathlib import Path

import pytest

from nf_core.modules.modules_json import ModulesJson
from nf_core.subworkflows.install import SubworkflowInstall

from ..test_subworkflows import TestSubworkflows
from ..utils import (
    CROSS_ORGANIZATION_URL,
    GITLAB_BRANCH_TEST_BRANCH,
    GITLAB_REPO,
    GITLAB_SUBWORKFLOWS_BRANCH,
    GITLAB_SUBWORKFLOWS_ORG_PATH_BRANCH,
    GITLAB_URL,
    with_temporary_folder,
)


class TestSubworkflowsInstall(TestSubworkflows):
    def test_subworkflows_install_bam_sort_stats_samtools(self):
        """Test installing a subworkflow - bam_sort_stats_samtools"""
        assert self.subworkflow_install.install("bam_sort_stats_samtools") is not False
        subworkflow_path = Path(
            self.subworkflow_install.directory, "subworkflows", "nf-core", "bam_sort_stats_samtools"
        )
        sub_subworkflow_path = Path(self.subworkflow_install.directory, "subworkflows", "nf-core", "bam_stats_samtools")
        samtools_index_path = Path(self.subworkflow_install.directory, "modules", "nf-core", "samtools", "index")
        samtools_sort_path = Path(self.subworkflow_install.directory, "modules", "nf-core", "samtools", "sort")
        samtools_stats_path = Path(self.subworkflow_install.directory, "modules", "nf-core", "samtools", "stats")
        samtools_idxstats_path = Path(self.subworkflow_install.directory, "modules", "nf-core", "samtools", "idxstats")
        samtools_flagstat_path = Path(self.subworkflow_install.directory, "modules", "nf-core", "samtools", "flagstat")
        assert subworkflow_path.exists()
        assert sub_subworkflow_path.exists()
        assert samtools_index_path.exists()
        assert samtools_sort_path.exists()
        assert samtools_stats_path.exists()
        assert samtools_idxstats_path.exists()
        assert samtools_flagstat_path.exists()

    def test_subworkflow_install_nopipeline(self):
        """Test installing a subworkflow - no pipeline given"""
        assert self.subworkflow_install.directory is not None
        self.subworkflow_install.directory = Path("non_existent_dir")
        assert self.subworkflow_install.install("bam_stats_samtools") is False

    @with_temporary_folder
    def test_subworkflows_install_emptypipeline(self, tmpdir):
        """Test installing a subworkflow - empty dir given"""

        Path(tmpdir, "nf-core-pipe").mkdir(exist_ok=True)
        self.subworkflow_install.directory = Path(tmpdir, "nf-core-pipe")
        with pytest.raises(UserWarning) as excinfo:
            self.subworkflow_install.install("bam_stats_samtools")
        assert "Could not find a 'main.nf' or 'nextflow.config' file" in str(excinfo.value)

    def test_subworkflows_install_nosubworkflow(self):
        """Test installing a subworkflow - unrecognised subworkflow given"""
        with pytest.raises(ValueError) as excinfo:
            self.subworkflow_install.install("foo")
        assert excinfo.typename == "ValueError"
        assert "Subworkflow 'foo' not found in available subworkflows" in self.caplog.text

    def test_subworkflows_install_bam_sort_stats_samtools_twice(self):
        """Test installing a subworkflow - bam_sort_stats_samtools already there"""
        self.subworkflow_install.install("bam_sort_stats_samtools")
        assert self.subworkflow_install.install("bam_sort_stats_samtools") is False

    def test_subworkflows_install_from_gitlab(self):
        """Test installing a subworkflow from GitLab"""
        assert self.subworkflow_install_gitlab.install("bam_stats_samtools") is True
        # Verify that the branch entry was added correctly
        modules_json = ModulesJson(self.pipeline_dir)
        assert (
            modules_json.get_component_branch(self.component_type, "bam_stats_samtools", GITLAB_URL, GITLAB_REPO)
            == GITLAB_SUBWORKFLOWS_BRANCH
        )

    def test_subworkflows_install_different_branch_fail(self):
        """Test installing a subworkflow from a different branch"""
        install_obj = SubworkflowInstall(self.pipeline_dir, remote_url=GITLAB_URL, branch=GITLAB_BRANCH_TEST_BRANCH)
        # The bam_stats_samtools subworkflow does not exists in the branch-test branch
        with pytest.raises(Exception) as excinfo:
            install_obj.install("bam_stats_samtools")
            assert "Subworkflow 'bam_stats_samtools' not found in available subworkflows" in str(excinfo.value)

    def test_subworkflows_install_across_organizations(self):
        """Test installing a subworkflow with modules from different organizations"""
        # The fastq_trim_fastp_fastqc subworkflow contains modules from different organizations
        self.subworkflow_install_cross_org.install("fastq_trim_fastp_fastqc")
        # Verify that the installed_by entry was added correctly
        modules_json = ModulesJson(self.pipeline_dir)
        mod_json = modules_json.get_modules_json()
        assert mod_json["repos"][CROSS_ORGANIZATION_URL]["modules"]["nf-core-test"]["fastqc"]["installed_by"] == [
            "fastq_trim_fastp_fastqc"
        ]

    def test_subworkflows_install_across_organizations_only_nfcore(self):
        """Test installing a subworkflow from a different organization but only with modules from nf-core"""
        # The wget_subwf subworkflow contains nf-core/modules/wget (and only that). It used to not be possible to install it twice
        # because of some org/nf-core confusion (https://github.com/nf-core/tools/issues/3876)
        # Note that we use two SubworkflowInstall to avoid cross-contamination
        for swfi in (self.subworkflow_install_cross_org, self.subworkflow_install_cross_org_again):
            swfi.install("wget_subwf")
            # Verify that the installed_by entry was added correctly
            modules_json = ModulesJson(self.pipeline_dir)
            mod_json = modules_json.get_modules_json()
            assert mod_json["repos"][CROSS_ORGANIZATION_URL]["subworkflows"]["nf-core-test"]["wget_subwf"][
                "installed_by"
            ] == ["subworkflows"]
            # This assertion used to fail on the second attempt
            assert (
                "wget_subwf"
                not in mod_json["repos"]["https://github.com/nf-core/modules.git"]["subworkflows"]["nf-core"]
            )

    def test_subworkflow_install_with_same_module(self):
        """Test installing a subworkflow with a module from a different organization that is already installed from another org"""
        # The fastq_trim_fastp_fastqc subworkflow contains the cross-org fastqc module, not the nf-core one
        self.subworkflow_install_cross_org.install("fastq_trim_fastp_fastqc")
        # Verify that the installed_by entry was added correctly
        modules_json = ModulesJson(self.pipeline_dir)
        mod_json = modules_json.get_modules_json()

        assert mod_json["repos"]["https://github.com/nf-core/modules.git"]["modules"]["nf-core"]["fastqc"][
            "installed_by"
        ] == ["modules"]

        assert mod_json["repos"][CROSS_ORGANIZATION_URL]["modules"]["nf-core-test"]["fastqc"]["installed_by"] == [
            "fastq_trim_fastp_fastqc"
        ]

    def test_subworkflows_install_tracking(self):
        """Test installing a subworkflow and finding the correct entries in installed_by section of modules.json"""
        assert self.subworkflow_install.install("bam_sort_stats_samtools")

        # Verify that the installed_by entry was added correctly
        modules_json = ModulesJson(self.pipeline_dir)
        mod_json = modules_json.get_modules_json()
        assert mod_json["repos"]["https://github.com/nf-core/modules.git"]["subworkflows"]["nf-core"][
            "bam_sort_stats_samtools"
        ]["installed_by"] == ["subworkflows"]
        assert mod_json["repos"]["https://github.com/nf-core/modules.git"]["subworkflows"]["nf-core"][
            "bam_stats_samtools"
        ]["installed_by"] == ["bam_sort_stats_samtools"]
        assert mod_json["repos"]["https://github.com/nf-core/modules.git"]["modules"]["nf-core"]["samtools/stats"][
            "installed_by"
        ] == ["bam_stats_samtools"]
        assert mod_json["repos"]["https://github.com/nf-core/modules.git"]["modules"]["nf-core"]["samtools/sort"][
            "installed_by"
        ] == ["bam_sort_stats_samtools"]

        # Clean directory
        self.subworkflow_remove.remove("bam_sort_stats_samtools")

    def test_subworkflows_install_tracking_added_already_installed(self):
        """Test installing a subworkflow and finding the correct entries in installed_by section of modules.json"""
        self.subworkflow_install.install("bam_sort_stats_samtools")
        self.subworkflow_install.install("bam_stats_samtools")

        # Verify that the installed_by entry was added correctly
        modules_json = ModulesJson(self.pipeline_dir)
        mod_json = modules_json.get_modules_json()
        assert mod_json["repos"]["https://github.com/nf-core/modules.git"]["subworkflows"]["nf-core"][
            "bam_sort_stats_samtools"
        ]["installed_by"] == ["subworkflows"]
        assert sorted(
            mod_json["repos"]["https://github.com/nf-core/modules.git"]["subworkflows"]["nf-core"][
                "bam_stats_samtools"
            ]["installed_by"]
        ) == sorted(["bam_sort_stats_samtools", "subworkflows"])

        # Clean directory
        self.subworkflow_remove.remove("bam_sort_stats_samtools")
        self.subworkflow_remove.remove("bam_stats_samtools")

    def test_subworkflows_install_tracking_added_super_subworkflow(self):
        """Test installing a subworkflow and finding the correct entries in installed_by section of modules.json"""
        self.subworkflow_install.install("bam_stats_samtools")
        self.subworkflow_install.install("bam_sort_stats_samtools")

        # Verify that the installed_by entry was added correctly
        modules_json = ModulesJson(self.pipeline_dir)
        mod_json = modules_json.get_modules_json()
        assert mod_json["repos"]["https://github.com/nf-core/modules.git"]["subworkflows"]["nf-core"][
            "bam_sort_stats_samtools"
        ]["installed_by"] == ["subworkflows"]
        assert sorted(
            mod_json["repos"]["https://github.com/nf-core/modules.git"]["subworkflows"]["nf-core"][
                "bam_stats_samtools"
            ]["installed_by"]
        ) == sorted(["subworkflows", "bam_sort_stats_samtools"])

    def test_subworkflows_install_alternate_remote(self):
        """Test installing a module from a different remote with the same organization path"""
        install_obj = SubworkflowInstall(
            self.pipeline_dir, remote_url=GITLAB_URL, branch=GITLAB_SUBWORKFLOWS_ORG_PATH_BRANCH
        )
        # Install a subworkflow from GitLab which is also installed from GitHub with the same org_path
        with pytest.raises(Exception) as excinfo:
            install_obj.install("fastqc")
            assert "Could not find a 'main.nf' or 'nextflow.config' file" in str(excinfo.value)
