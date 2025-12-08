from pathlib import Path

from nf_core.modules.modules_json import ModulesJson
from tests.utils import CROSS_ORGANIZATION_URL

from ..test_subworkflows import TestSubworkflows


class TestSubworkflowsRemove(TestSubworkflows):
    def test_subworkflows_remove_uninstalled_subworkflow(self):
        """Test removing subworkflow without installing it"""
        assert self.subworkflow_remove.remove("bam_sort_stats_samtools") is False

    def test_subworkflows_remove_subworkflow(self):
        """Test removing subworkflow and all it's dependencies after installing it"""
        self.subworkflow_install.install("bam_sort_stats_samtools")

        subworkflow_path = Path(self.subworkflow_install.directory, "subworkflows", "nf-core")
        bam_sort_stats_samtools_path = Path(subworkflow_path, "bam_sort_stats_samtools")
        bam_stats_samtools_path = Path(subworkflow_path, "bam_stats_samtools")
        samtools_index_path = Path(self.subworkflow_install.directory, "modules", "nf-core", "samtools", "index")
        ModulesJson(self.pipeline_dir)
        mod_json_before = ModulesJson(self.pipeline_dir).get_modules_json()
        assert self.subworkflow_remove.remove("bam_sort_stats_samtools")
        mod_json_after = ModulesJson(self.pipeline_dir).get_modules_json()
        assert Path.exists(bam_sort_stats_samtools_path) is False
        assert Path.exists(bam_stats_samtools_path) is False
        assert Path.exists(samtools_index_path) is False
        assert mod_json_before != mod_json_after
        # assert subworkflows key is removed from modules.json
        assert (
            "bam_sort_stats_samtools"
            not in mod_json_after["repos"]["https://github.com/nf-core/modules.git"]["subworkflows"].keys()
        )
        assert (
            "samtools/index" not in mod_json_after["repos"]["https://github.com/nf-core/modules.git"]["modules"].keys()
        )

    def test_subworkflows_remove_subworkflow_keep_installed_module(self):
        """Test removing subworkflow and all it's dependencies after installing it, except for a separately installed module"""
        self.subworkflow_install.install("bam_sort_stats_samtools")
        self.mods_install.install("samtools/index")

        subworkflow_path = Path(self.subworkflow_install.directory, "subworkflows", "nf-core")
        bam_sort_stats_samtools_path = Path(subworkflow_path, "bam_sort_stats_samtools")
        bam_stats_samtools_path = Path(subworkflow_path, "bam_stats_samtools")
        samtools_index_path = Path(self.subworkflow_install.directory, "modules", "nf-core", "samtools", "index")

        mod_json_before = ModulesJson(self.pipeline_dir).get_modules_json()
        assert self.subworkflow_remove.remove("bam_sort_stats_samtools")
        mod_json_after = ModulesJson(self.pipeline_dir).get_modules_json()

        assert Path.exists(bam_sort_stats_samtools_path) is False
        assert Path.exists(bam_stats_samtools_path) is False
        assert Path.exists(samtools_index_path) is True
        assert mod_json_before != mod_json_after
        # assert subworkflows key is removed from modules.json
        assert (
            "bam_sort_stats_samtools"
            not in mod_json_after["repos"]["https://github.com/nf-core/modules.git"]["subworkflows"].keys()
        )
        assert (
            "samtools/index"
            in mod_json_after["repos"]["https://github.com/nf-core/modules.git"]["modules"]["nf-core"].keys()
        )

    def test_subworkflows_remove_one_of_two_subworkflow(self):
        """Test removing subworkflow and all it's dependencies after installing it"""
        self.subworkflow_install.install("bam_sort_stats_samtools")
        self.subworkflow_install.install("bam_stats_samtools")
        subworkflow_path = Path(self.subworkflow_install.directory, "subworkflows", "nf-core")
        bam_sort_stats_samtools_path = Path(subworkflow_path, "bam_sort_stats_samtools")
        bam_stats_samtools_path = Path(subworkflow_path, "bam_stats_samtools")
        samtools_index_path = Path(self.subworkflow_install.directory, "modules", "nf-core", "samtools", "index")
        samtools_stats_path = Path(self.subworkflow_install.directory, "modules", "nf-core", "samtools", "stats")

        assert self.subworkflow_remove.remove("bam_sort_stats_samtools")

        assert Path.exists(subworkflow_path) is True
        assert Path.exists(bam_sort_stats_samtools_path) is False
        assert Path.exists(bam_stats_samtools_path) is True
        assert Path.exists(samtools_index_path) is False
        assert Path.exists(samtools_stats_path) is True
        self.subworkflow_remove.remove("bam_stats_samtools")

    def test_subworkflows_remove_included_subworkflow(self):
        """Test removing subworkflow which is installed by another subworkflow and all it's dependencies."""
        self.subworkflow_install.install("bam_sort_stats_samtools")
        subworkflow_path = Path(self.subworkflow_install.directory, "subworkflows", "nf-core")
        bam_sort_stats_samtools_path = Path(subworkflow_path, "bam_sort_stats_samtools")
        bam_stats_samtools_path = Path(subworkflow_path, "bam_stats_samtools")
        samtools_index_path = Path(self.subworkflow_install.directory, "modules", "nf-core", "samtools", "index")
        samtools_stats_path = Path(self.subworkflow_install.directory, "modules", "nf-core", "samtools", "stats")

        assert self.subworkflow_remove.remove("bam_stats_samtools") is False

        assert Path.exists(subworkflow_path) is True
        assert Path.exists(bam_sort_stats_samtools_path) is True
        assert Path.exists(bam_stats_samtools_path) is True
        assert Path.exists(samtools_index_path) is True
        assert Path.exists(samtools_stats_path) is True
        self.subworkflow_remove.remove("bam_sort_stats_samtools")

    def test_subworkflows_remove_subworkflow_keep_installed_cross_org_module(self):
        """Test removing subworkflow and all it's dependencies after installing it, except for a separately installed module from another organisation"""
        self.subworkflow_install_cross_org.install("fastq_trim_fastp_fastqc")
        self.mods_install.install("fastqc")

        subworkflow_path = Path(self.subworkflow_install.directory, "subworkflows", "nf-core-test")
        fastq_trim_fastp_fastqc_path = Path(subworkflow_path, "fastq_trim_fastp_fastqc")
        fastqc_path = Path(self.subworkflow_install.directory, "modules", "nf-core-test", "fastqc")
        nfcore_fastqc_path = Path(self.subworkflow_install.directory, "modules", "nf-core", "fastqc")

        mod_json_before = ModulesJson(self.pipeline_dir).get_modules_json()
        assert self.subworkflow_remove_cross_org.remove("fastq_trim_fastp_fastqc")
        mod_json_after = ModulesJson(self.pipeline_dir).get_modules_json()

        assert Path.exists(fastq_trim_fastp_fastqc_path) is False
        assert Path.exists(fastqc_path) is False
        assert Path.exists(nfcore_fastqc_path) is True
        assert mod_json_before != mod_json_after
        # assert subworkflows key is removed from modules.json
        assert CROSS_ORGANIZATION_URL not in mod_json_after["repos"].keys()
        assert (
            "fastqc" in mod_json_after["repos"]["https://github.com/nf-core/modules.git"]["modules"]["nf-core"].keys()
        )
        assert (
            "fastp"
            not in mod_json_after["repos"]["https://github.com/nf-core/modules.git"]["modules"]["nf-core"].keys()
        )
