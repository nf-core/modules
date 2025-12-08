import json
import shutil
import subprocess
from pathlib import Path

import nf_core.subworkflows

from ..test_subworkflows import TestSubworkflows
from ..utils import GITLAB_SUBWORKFLOWS_BRANCH, GITLAB_URL


class TestSubworkflowsLint(TestSubworkflows):
    def test_subworkflows_lint(self):
        """Test linting the fastq_align_bowtie2 subworkflow"""
        self.subworkflow_install.install("fastq_align_bowtie2")
        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.pipeline_dir)
        subworkflow_lint.lint(print_results=False, subworkflow="fastq_align_bowtie2")
        assert len(subworkflow_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in subworkflow_lint.failed]}"
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) >= 0

    def test_subworkflows_lint_empty(self):
        """Test linting a pipeline with no subworkflows installed"""
        self.subworkflow_remove.remove("utils_nextflow_pipeline", force=True)
        self.subworkflow_remove.remove("utils_nfcore_pipeline", force=True)
        self.subworkflow_remove.remove("utils_nfschema_plugin", force=True)
        nf_core.subworkflows.SubworkflowLint(directory=self.pipeline_dir)
        assert "No subworkflows from https://github.com/nf-core/modules.git installed in pipeline" in self.caplog.text

    def test_subworkflows_lint_new_subworkflow(self):
        """lint a new subworkflow"""
        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.nfcore_modules)
        subworkflow_lint.lint(print_results=True, all_subworkflows=True)
        assert len(subworkflow_lint.failed) == 0
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) >= 0

    def test_subworkflows_lint_no_gitlab(self):
        """Test linting a pipeline with no subworkflows installed"""
        nf_core.subworkflows.SubworkflowLint(directory=self.pipeline_dir, remote_url=GITLAB_URL)
        assert f"No subworkflows from {GITLAB_URL} installed in pipeline" in self.caplog.text

    def test_subworkflows_lint_gitlab_subworkflows(self):
        """Lint subworkflows from a different remote"""
        self.subworkflow_install_gitlab.install("bam_stats_samtools")
        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(
            directory=self.pipeline_dir, remote_url=GITLAB_URL, branch=GITLAB_SUBWORKFLOWS_BRANCH
        )
        subworkflow_lint.lint(print_results=False, all_subworkflows=True)
        assert len(subworkflow_lint.failed) == 0
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) >= 0

    def test_subworkflows_lint_multiple_remotes(self):
        """Lint subworkflows from a different remote"""
        self.subworkflow_install_gitlab.install("bam_stats_samtools")
        self.subworkflow_install.install("fastq_align_bowtie2")
        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(
            directory=self.pipeline_dir, remote_url=GITLAB_URL, branch=GITLAB_SUBWORKFLOWS_BRANCH
        )
        subworkflow_lint.lint(print_results=False, all_subworkflows=True)
        assert len(subworkflow_lint.failed) == 0
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) >= 0

    def test_subworkflows_lint_update_meta_yml(self):
        """update the meta.yml of a subworkflow"""
        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.nfcore_modules, fix=True)
        subworkflow_lint.lint(print_results=False, subworkflow="test_subworkflow")
        assert len(subworkflow_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in subworkflow_lint.failed]}"
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) >= 0

    def test_subworkflows_lint_snapshot_file(self):
        """Test linting a subworkflow with a snapshot file"""
        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.nfcore_modules)
        subworkflow_lint.lint(print_results=False, subworkflow="test_subworkflow")
        assert len(subworkflow_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in subworkflow_lint.failed]}"
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) >= 0

    def test_subworkflows_lint_snapshot_file_missing_fail(self):
        """Test linting a subworkflow with a snapshot file missing, which should fail"""
        Path(
            self.nfcore_modules,
            "subworkflows",
            "nf-core",
            "test_subworkflow",
            "tests",
            "main.nf.test.snap",
        ).unlink()
        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.nfcore_modules)
        subworkflow_lint.lint(print_results=False, subworkflow="test_subworkflow")
        Path(
            self.nfcore_modules,
            "subworkflows",
            "nf-core",
            "test_subworkflow",
            "tests",
            "main.nf.test.snap",
        ).touch()
        assert len(subworkflow_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in subworkflow_lint.failed]}"
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) >= 0

    def test_subworkflows_lint_snapshot_file_not_needed(self):
        """Test linting a subworkflow which doesn't need a snapshot file by removing the snapshot keyword in the main.nf.test file"""
        with open(
            Path(
                self.nfcore_modules,
                "subworkflows",
                "nf-core",
                "test_subworkflow",
                "tests",
                "main.nf.test",
            )
        ) as fh:
            content = fh.read()
            new_content = content.replace("snapshot(", "snap (")
        with open(
            Path(
                self.nfcore_modules,
                "subworkflows",
                "nf-core",
                "test_subworkflow",
                "tests",
                "main.nf.test",
            ),
            "w",
        ) as fh:
            fh.write(new_content)

        Path(
            self.nfcore_modules,
            "subworkflows",
            "nf-core",
            "test_subworkflow",
            "tests",
            "main.nf.test.snap",
        ).unlink()
        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.nfcore_modules)
        subworkflow_lint.lint(print_results=False, subworkflow="test_subworkflow")
        Path(
            self.nfcore_modules,
            "subworkflows",
            "nf-core",
            "test_subworkflow",
            "tests",
            "main.nf.test.snap",
        ).touch()
        assert len(subworkflow_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in subworkflow_lint.failed]}"
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) >= 0

    def test_subworkflows_lint_less_than_two_modules_warning(self):
        """Test linting a subworkflow with less than two modules"""
        self.subworkflow_install.install("bam_stats_samtools")
        # Remove two modules
        with open(
            Path(
                self.pipeline_dir,
                "subworkflows",
                "nf-core",
                "bam_stats_samtools",
                "main.nf",
            )
        ) as fh:
            content = fh.read()
            new_content = content.replace(
                "include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'",
                "",
            )
            new_content = new_content.replace(
                "include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'",
                "",
            )
        with open(
            Path(
                self.pipeline_dir,
                "subworkflows",
                "nf-core",
                "bam_stats_samtools",
                "main.nf",
            ),
            "w",
        ) as fh:
            fh.write(new_content)
        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.pipeline_dir)
        subworkflow_lint.lint(print_results=False, subworkflow="bam_stats_samtools")
        assert len(subworkflow_lint.failed) >= 0, f"Linting failed with {[x.__dict__ for x in subworkflow_lint.failed]}"
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) > 0
        assert subworkflow_lint.warned[0].lint_test == "main_nf_include"
        # cleanup
        self.subworkflow_remove.remove("bam_stats_samtools", force=True)

    def test_subworkflows_lint_include_multiple_alias(self):
        """Test linting a subworkflow with multiple include methods"""
        self.subworkflow_install.install("bam_stats_samtools")
        with open(
            Path(
                self.pipeline_dir,
                "subworkflows",
                "nf-core",
                "bam_stats_samtools",
                "main.nf",
            )
        ) as fh:
            content = fh.read()
            new_content = content.replace("SAMTOOLS_STATS", "SAMTOOLS_STATS_1")
            new_content = new_content.replace(
                "include { SAMTOOLS_STATS_1 ",
                "include { SAMTOOLS_STATS as SAMTOOLS_STATS_1; SAMTOOLS_STATS as SAMTOOLS_STATS_2 ",
            )
        with open(
            Path(
                self.pipeline_dir,
                "subworkflows",
                "nf-core",
                "bam_stats_samtools",
                "main.nf",
            ),
            "w",
        ) as fh:
            fh.write(new_content)

        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.pipeline_dir)
        subworkflow_lint.lint(print_results=False, subworkflow="bam_stats_samtools")
        assert len(subworkflow_lint.failed) >= 0, f"Linting failed with {[x.__dict__ for x in subworkflow_lint.failed]}"
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) == 2
        assert any(
            [
                x.message == "Included component 'SAMTOOLS_STATS_1' versions are added in main.nf"
                for x in subworkflow_lint.passed
            ]
        )
        assert any(
            [x.message == "Included component 'SAMTOOLS_STATS_1' used in main.nf" for x in subworkflow_lint.passed]
        )
        assert any(
            [x.message == "Included component 'SAMTOOLS_STATS_2' not used in main.nf" for x in subworkflow_lint.warned]
        )

        # cleanup
        self.subworkflow_remove.remove("bam_stats_samtools", force=True)

    def test_subworkflows_lint_capitalization_fail(self):
        """Test linting a subworkflow with a capitalization fail"""
        self.subworkflow_install.install("bam_stats_samtools")
        # change workflow name to lowercase
        with open(
            Path(
                self.pipeline_dir,
                "subworkflows",
                "nf-core",
                "bam_stats_samtools",
                "main.nf",
            )
        ) as fh:
            content = fh.read()
            new_content = content.replace("workflow BAM_STATS_SAMTOOLS {", "workflow bam_stats_samtools {")
        with open(
            Path(
                self.pipeline_dir,
                "subworkflows",
                "nf-core",
                "bam_stats_samtools",
                "main.nf",
            ),
            "w",
        ) as fh:
            fh.write(new_content)
        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.pipeline_dir)
        subworkflow_lint.lint(print_results=False, subworkflow="bam_stats_samtools")
        assert len(subworkflow_lint.failed) >= 1, f"Linting failed with {[x.__dict__ for x in subworkflow_lint.failed]}"
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) >= 0
        assert any([x.lint_test == "workflow_capitals" for x in subworkflow_lint.failed])

        # cleanup
        self.subworkflow_remove.remove("bam_stats_samtools", force=True)

    def test_subworkflows_absent_version(self):
        """Test linting a nf-test subworkflow if the versions is absent in the snapshot file `"""
        snap_file = Path(
            self.nfcore_modules,
            "subworkflows",
            "nf-core",
            "test_subworkflow",
            "tests",
            "main.nf.test.snap",
        )
        with open(snap_file) as fh:
            content = fh.read()
            new_content = content.replace("versions", "foo")
        with open(snap_file, "w") as fh:
            fh.write(new_content)

        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.nfcore_modules)
        subworkflow_lint.lint(print_results=False, subworkflow="test_subworkflow")
        assert len(subworkflow_lint.failed) == 0
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) >= 0, f"Linting warned with {[x.__dict__ for x in subworkflow_lint.warned]}"
        assert any([x.lint_test == "test_snap_versions" for x in subworkflow_lint.warned])

        # cleanup
        with open(snap_file, "w") as fh:
            fh.write(content)

    def test_subworkflows_missing_test_dir(self):
        """Test linting a nf-test subworkflow if the tests directory is missing"""
        test_dir = Path(self.nfcore_modules, "subworkflows", "nf-core", "test_subworkflow", "tests")
        test_dir_copy = shutil.copytree(test_dir, test_dir.parent / "tests_copy")
        shutil.rmtree(test_dir)

        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(self.nfcore_modules)
        subworkflow_lint.lint(print_results=False, subworkflow="test_subworkflow")
        assert len(subworkflow_lint.failed) == 1
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) >= 0, f"Linting warned with {[x.__dict__ for x in subworkflow_lint.warned]}"
        assert any([x.lint_test == "test_dir_exists" for x in subworkflow_lint.failed])

        # cleanup
        shutil.copytree(test_dir_copy, test_dir)

    # There are many steps before the actual main_nf linting where we rely on the main_nf file to exist, so this test is not possible for now
    # def test_subworkflows_missing_main_nf(self):
    #     """Test linting a nf-test subworkflow if the main.nf file is missing"""

    #     subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.nfcore_modules)
    #     main_nf = Path(self.nfcore_modules, "subworkflows", "nf-core", "test_subworkflow", "main.nf")
    #     main_nf_copy = shutil.copy(main_nf, main_nf.parent / "main_nf_copy")
    #     main_nf.unlink()
    #     subworkflow_lint.lint(print_results=False, subworkflow="test_subworkflow")
    #     assert len(subworkflow_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in subworkflow_lint.failed]}"
    #     assert len(subworkflow_lint.passed) > 0
    #     assert len(subworkflow_lint.warned) >= 0
    #     assert subworkflow_lint.failed[0].lint_test == "main_nf_exists"

    #     # cleanup
    #     shutil.copy(main_nf_copy, main_nf)
    #     shutil.rmtree(Path(self.nfcore_modules, "subworkflows", "nf-core", "test_subworkflow_backup"))

    def test_subworkflows_empty_file_in_snapshot(self):
        """Test linting a nf-test subworkflow with an empty file sha sum in the test snapshot, which should make it fail (if it is not a stub)"""
        snap_file = Path(
            self.nfcore_modules,
            "subworkflows",
            "nf-core",
            "test_subworkflow",
            "tests",
            "main.nf.test.snap",
        )
        snap = json.load(snap_file.open())
        content = snap_file.read_text()
        snap["my test"]["content"][0]["0"] = "test:md5,d41d8cd98f00b204e9800998ecf8427e"

        with open(snap_file, "w") as fh:
            json.dump(snap, fh)

        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.nfcore_modules)
        subworkflow_lint.lint(print_results=False, subworkflow="test_subworkflow")
        assert len(subworkflow_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in subworkflow_lint.failed]}"
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) >= 0
        assert subworkflow_lint.failed[0].lint_test == "test_snap_md5sum"

        # reset the file
        with open(snap_file, "w") as fh:
            fh.write(content)

    def test_subworkflows_empty_file_in_stub_snapshot(self):
        """Test linting a nf-test subworkflow with an empty file sha sum in the stub test snapshot, which should make it not fail"""
        snap_file = Path(
            self.nfcore_modules,
            "subworkflows",
            "nf-core",
            "test_subworkflow",
            "tests",
            "main.nf.test.snap",
        )
        snap = json.load(snap_file.open())
        content = snap_file.read_text()
        snap["my_test_stub"] = {"content": [{"0": "test:md5,d41d8cd98f00b204e9800998ecf8427e", "versions": {}}]}

        with open(snap_file, "w") as fh:
            json.dump(snap, fh)

        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.nfcore_modules)
        subworkflow_lint.lint(print_results=False, subworkflow="test_subworkflow")
        assert len(subworkflow_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in subworkflow_lint.failed]}"
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) >= 0
        assert any(x.lint_test == "test_snap_md5sum" for x in subworkflow_lint.passed)

        # reset the file
        with open(snap_file, "w") as fh:
            fh.write(content)

    def test_subworkflows_lint_local(self):
        assert self.subworkflow_install.install("fastq_align_bowtie2")
        installed = Path(self.pipeline_dir, "subworkflows", "nf-core", "fastq_align_bowtie2")
        local = Path(self.pipeline_dir, "subworkflows", "local", "fastq_align_bowtie2")
        shutil.move(installed, local)
        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.pipeline_dir)
        subworkflow_lint.lint(print_results=False, local=True)
        assert len(subworkflow_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in subworkflow_lint.failed]}"
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) >= 0

    def test_subworkflows_lint_local_missing_files(self):
        assert self.subworkflow_install.install("fastq_align_bowtie2")
        installed = Path(self.pipeline_dir, "subworkflows", "nf-core", "fastq_align_bowtie2")
        local = Path(self.pipeline_dir, "subworkflows", "local", "fastq_align_bowtie2")
        shutil.move(installed, local)
        Path(self.pipeline_dir, "subworkflows", "local", "fastq_align_bowtie2", "meta.yml").unlink()
        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.pipeline_dir)
        subworkflow_lint.lint(print_results=False, local=True)
        assert len(subworkflow_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in subworkflow_lint.failed]}"
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) >= 0
        warnings = [x.message for x in subworkflow_lint.warned]
        assert "Subworkflow `meta.yml` does not exist" in warnings

    def test_subworkflows_lint_local_old_format(self):
        assert self.subworkflow_install.install("fastq_align_bowtie2")
        installed = Path(self.pipeline_dir, "subworkflows", "nf-core", "fastq_align_bowtie2", "main.nf")
        Path(self.pipeline_dir, "subworkflows", "local").mkdir(exist_ok=True)
        local = Path(self.pipeline_dir, "subworkflows", "local", "fastq_align_bowtie2.nf")
        shutil.copy(installed, local)
        self.subworkflow_remove.remove("fastq_align_bowtie2", force=True)
        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.pipeline_dir)
        subworkflow_lint.lint(print_results=False, local=True)
        assert len(subworkflow_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in subworkflow_lint.failed]}"
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) >= 0


class TestSubworkflowsLintPatch(TestSubworkflows):
    def setUp(self) -> None:
        super().setUp()

        # Install the subworkflow bam_sort_stats_samtools
        self.subworkflow_install.install("bam_sort_stats_samtools")

        # Modify the subworkflow by inserting a new input channel
        new_line = "    ch_dummy // channel: [ path ]\n"

        subworkflow_path = Path(self.pipeline_dir, "subworkflows", "nf-core", "bam_sort_stats_samtools", "main.nf")

        with open(subworkflow_path) as fh:
            lines = fh.readlines()
        for line_index in range(len(lines)):
            if "take:" in lines[line_index]:
                lines.insert(line_index + 1, new_line)
        with open(subworkflow_path, "w") as fh:
            fh.writelines(lines)

        # Create a patch file
        self.patch_obj = nf_core.subworkflows.SubworkflowPatch(self.pipeline_dir)
        self.patch_obj.patch("bam_sort_stats_samtools")

    def test_lint_clean_patch(self):
        """Test linting a patched subworkflow"""

        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.pipeline_dir)
        subworkflow_lint.lint(print_results=False, subworkflow="bam_sort_stats_samtools")

        assert len(subworkflow_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in subworkflow_lint.failed]}"
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) == 0, f"Linting warned with {[x.__dict__ for x in subworkflow_lint.warned]}"

    def test_lint_broken_patch(self):
        """Test linting a patched subworkflow when the patch is broken"""

        # Now modify the diff
        diff_file = Path(
            self.pipeline_dir, "subworkflows", "nf-core", "bam_sort_stats_samtools", "bam_sort_stats_samtools.diff"
        )
        subprocess.check_call(["sed", "-i''", "s/...$//", str(diff_file)])

        subworkflow_lint = nf_core.subworkflows.SubworkflowLint(directory=self.pipeline_dir)
        subworkflow_lint.lint(print_results=False, subworkflow="bam_sort_stats_samtools")

        assert len(subworkflow_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in subworkflow_lint.failed]}"
        errors = [x.message for x in subworkflow_lint.failed]
        assert "Subworkflow patch cannot be cleanly applied" in errors
        assert len(subworkflow_lint.passed) > 0
        assert len(subworkflow_lint.warned) == 0, f"Linting warned with {[x.__dict__ for x in subworkflow_lint.warned]}"
