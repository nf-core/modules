process SAMTOOLS_GETRG {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c2f02ceb8cbe5aa37930dcbc0de6f1d07dfb292178c7e780c44ab2f556bf829f/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--8a3ddf262165ba23'}"

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
