process BAM2FASTX_BAM2FASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bam2fastx:1.3.1--hf05d43a_1':
        'quay.io/biocontainers/bam2fastx:1.3.1--hf05d43a_1' }"

    input:
    tuple val(meta), path(bam), path(index)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    tuple val("${task.process}"), val('bam2fastx'), eval("bam2fastq --version 2>&1) | sed 's/^.*bam2fastq //'"), emit: versions_bam2fastx, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/pbtk/bam2fastq

Reason:
This module is no longer fit for purpose because bam2fastx has been deprecated by PacificBiosciences

"""
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    assert false: deprecation_message
    """
    bam2fastq \\
        $args \\
        -o ${prefix} \\
        $bam
    """
}
