process SAMTOOLS_GETRG {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), file("readgroups.txt"), emit: readgroup
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/samtools/splitheader

Reason:
This module has been renamed to samtools/splitheader, which has the same functionality but
extends the outputs to include other types of SAM header.
"""
    assert false : deprecation_message
    """
    """

    stub:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/samtools/splitheader

Reason:
This module has been renamed to samtools/splitheader, which has the same functionality but
extends the outputs to include other types of SAM header.
"""
    assert false : deprecation_message
    """
    """
}
