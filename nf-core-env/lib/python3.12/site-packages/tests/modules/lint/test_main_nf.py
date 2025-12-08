import pytest

import nf_core.modules.lint
from nf_core.components.nfcore_component import NFCoreComponent
from nf_core.modules.lint.main_nf import check_container_link_line, check_process_labels

from ...test_modules import TestModules
from .test_lint_utils import MockModuleLint


@pytest.mark.parametrize(
    "content,passed,warned,failed",
    [
        # Valid process label
        ("label 'process_high'\ncpus 12", 1, 0, 0),
        # Non-alphanumeric characters in label
        ("label 'a:label:with:colons'\ncpus 12", 0, 2, 0),
        # Conflicting labels
        ("label 'process_high'\nlabel 'process_low'\ncpus 12", 0, 1, 0),
        # Duplicate labels
        ("label 'process_high'\nlabel 'process_high'\ncpus 12", 0, 2, 0),
        # Valid and non-standard labels
        ("label 'process_high'\nlabel 'process_extra_label'\ncpus 12", 1, 1, 0),
        # Non-standard label only
        ("label 'process_extra_label'\ncpus 12", 0, 2, 0),
        # Non-standard duplicates without quotes
        ("label process_extra_label\nlabel process_extra_label\ncpus 12", 0, 3, 0),
        # No label found
        ("cpus 12", 0, 1, 0),
    ],
)
def test_process_labels(content, passed, warned, failed):
    """Test process label validation"""
    mock_lint = MockModuleLint()
    check_process_labels(mock_lint, content.splitlines())

    assert len(mock_lint.passed) == passed
    assert len(mock_lint.warned) == warned
    assert len(mock_lint.failed) == failed


@pytest.mark.parametrize(
    "content,passed,warned,failed",
    [
        # Single-line container definition should pass
        ('container "quay.io/nf-core/gatk:4.4.0.0" //Biocontainers is missing a package', 2, 0, 0),
        # Multi-line container definition should pass
        (
            '''container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
                'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"''',
            6,
            0,
            0,
        ),
        # Space in container URL should fail
        (
            '''container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0 ':
                'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"''',
            5,
            0,
            1,
        ),
        # Incorrect quoting of container string should fail
        (
            '''container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0 ':
                "biocontainers/gatk4:4.4.0.0--py36hdfd78af_0" }"''',
            4,
            0,
            1,
        ),
        # Ternary with ? on next line (new Nextflow format) should pass
        (
            '''container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c262fc09eca59edb5a724080eeceb00fb06396f510aefb229c2d2c6897e63975/data'
        : 'community.wave.seqera.io/library/coreutils:9.5--ae99c88a9b28c264'}"''',
            6,
            0,
            0,
        ),
    ],
)
def test_container_links(content, passed, warned, failed):
    """Test container link validation"""
    mock_lint = MockModuleLint()

    for line in content.splitlines():
        if line.strip():
            check_container_link_line(mock_lint, line, registry="quay.io")

    assert len(mock_lint.passed) == passed
    assert len(mock_lint.warned) == warned
    assert len(mock_lint.failed) == failed


class TestMainNfLinting(TestModules):
    """
    Test main.nf linting functionality.

    This class tests various aspects of main.nf file linting including:
    - Process label validation and standards compliance
    - Container definition syntax and URL validation
    - Integration testing with alternative registries
    - General module linting workflow
    """

    def setUp(self):
        """Set up test fixtures by installing required modules"""
        super().setUp()
        # Install samtools/sort module for all tests in this class
        if not self.mods_install.install("samtools/sort"):
            self.skipTest("Could not install samtools/sort module")
        if not self.mods_install.install("bamstats/generalstats"):
            self.skipTest("Could not install samtools/sort module")

    def test_main_nf_lint_with_alternative_registry(self):
        """Test main.nf linting with alternative container registry"""
        # Test with alternative registry - should warn/fail when containers don't match the registry
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir, registry="public.ecr.aws")
        module_lint.lint(print_results=False, module="samtools/sort")

        # Alternative registry should produce warnings or failures for container mismatches
        # since samtools/sort module likely uses biocontainers/quay.io, not public.ecr.aws
        total_issues = len(module_lint.failed) + len(module_lint.warned)
        assert total_issues > 0, (
            "Expected warnings/failures when using alternative registry that doesn't match module containers"
        )

        # Test with default registry - should pass cleanly
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, module="samtools/sort")
        assert len(module_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) > 0

    def test_topics_and_emits_version_check(self):
        """Test that main_nf version emit and topics check works correctly"""

        # Lint a module known to have versions YAML in main.nf (for now)
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, module="samtools/sort")
        assert len(module_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.warned) == 2, (
            f"Linting warned with {[x.__dict__ for x in module_lint.warned]}, expected 2 warnings"
        )
        assert len(module_lint.passed) > 0

        # Lint a module known to have topics as output in main.nf
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.pipeline_dir)
        module_lint.lint(print_results=False, module="bamstats/generalstats")
        assert len(module_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.warned) == 0, (
            f"Linting warned with {[x.__dict__ for x in module_lint.warned]}, expected 1 warning"
        )
        assert len(module_lint.passed) > 0


def test_get_inputs_no_partial_keyword_match(tmp_path):
    """Test that input parsing doesn't match keywords within larger words like 'evaluate' or 'pathogen'"""
    main_nf_content = """
process TEST_PROCESS {
    input:
    val(meta)
    path(reads)
    tuple val(evaluate), path(pathogen)

    output:
    path("*.txt"), emit: results

    script:
    "echo test"
}
"""
    main_nf_path = tmp_path / "main.nf"
    main_nf_path.write_text(main_nf_content)

    component = NFCoreComponent(
        component_name="test",
        repo_url=None,
        component_dir=tmp_path,
        repo_type="modules",
        base_dir=tmp_path,
        component_type="modules",
        remote_component=False,
    )

    component.get_inputs_from_main_nf()

    # Should find 3 inputs: meta, reads, and the tuple (evaluate, pathogen)
    # The regex with \b should correctly identify 'val(evaluate)' and 'path(pathogen)' as valid inputs
    assert len(component.inputs) == 3, f"Expected 3 inputs, got {len(component.inputs)}: {component.inputs}"
    assert {"meta": {}} in component.inputs
    assert {"reads": {}} in component.inputs
    # The tuple should be captured as a list of two elements
    tuple_input = [{"evaluate": {}}, {"pathogen": {}}]
    assert tuple_input in component.inputs


def test_get_outputs_no_partial_keyword_match(tmp_path):
    """Test that output parsing doesn't match keywords within larger words like 'evaluate' or 'pathogen'"""
    main_nf_content = """
process TEST_PROCESS {
    input:
    val(meta)

    output:
    path("*.txt"), emit: results
    val(evaluate_result), emit: evaluation
    path(pathogen_data), emit: pathogens

    script:
    "echo test"
}
"""
    main_nf_path = tmp_path / "main.nf"
    main_nf_path.write_text(main_nf_content)

    component = NFCoreComponent(
        component_name="test",
        repo_url=None,
        component_dir=tmp_path,
        repo_type="modules",
        base_dir=tmp_path,
        component_type="modules",
        remote_component=False,
    )

    component.get_outputs_from_main_nf()

    # Should find 3 outputs with variable names containing 'val' and 'path' substrings
    # The regex with \b should correctly identify val(evaluate_result) and path(pathogen_data)
    assert len(component.outputs) == 3, f"Expected 3 outputs, got {len(component.outputs)}: {component.outputs}"
    assert "results" in component.outputs
    assert "evaluation" in component.outputs
    assert "pathogens" in component.outputs


def test_get_topics_no_partial_keyword_match(tmp_path):
    """Test that topic parsing doesn't match keywords within larger words like 'evaluate'"""
    main_nf_content = """
process TEST_PROCESS {
    input:
    val(meta)

    output:
    path("*.txt"), topic: results
    val(evaluate_result), topic: evaluation

    script:
    "echo test"
}
"""
    main_nf_path = tmp_path / "main.nf"
    main_nf_path.write_text(main_nf_content)

    component = NFCoreComponent(
        component_name="test",
        repo_url=None,
        component_dir=tmp_path,
        repo_type="modules",
        base_dir=tmp_path,
        component_type="modules",
        remote_component=False,
    )

    component.get_topics_from_main_nf()

    # Should find 2 topics with variable names containing 'val' substring
    # The regex with \b should correctly identify val(evaluate_result)
    assert len(component.topics) == 2, f"Expected 2 topics, got {len(component.topics)}: {component.topics}"
    assert "results" in component.topics
    assert "evaluation" in component.topics
