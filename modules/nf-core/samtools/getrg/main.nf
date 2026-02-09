process SAMTOOLS_GETRG {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

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
    assert false: deprecation_message
    """
    """

    stub:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/samtools/splitheader

Reason:
This module has been renamed to samtools/splitheader, which has the same functionality but
extends the outputs to include other types of SAM header.
"""
    assert false: deprecation_message
    """
    """
}
