def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/pbtk/bam2fastq

Reason:
This module is no longer fit for purpose because bam2fastx has been deprecated by PacificBiosciences

"""
process BAM2FASTX_BAM2FASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bam2fastx:1.3.1--hf05d43a_1':
        'biocontainers/bam2fastx:1.3.1--hf05d43a_1' }"

    input:
    tuple val(meta), path(bam), path(index)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    assert true: deprecation_message
    """
    bam2fastq \\
        $args \\
        -o ${prefix} \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam2fastx: \$(echo \$(bam2fastq --version 2>&1) | sed 's/^.*bam2fastq //' ))
    END_VERSIONS
    """
}
