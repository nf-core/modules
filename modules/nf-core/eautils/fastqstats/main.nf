process EAUTILS_FASTQSTATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ea-utils:1.1.2.779--h9dd4a16_0':
        'quay.io/biocontainers/ea-utils:1.1.2.779--h9dd4a16_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*_general.tsv") , emit: general_stats
    tuple val(meta), path("*_per_base.tsv"), emit: per_base_stats
    // tool does not have a version command, update this line on version bump
    tuple val("${task.process}"), val('fastq-stats'), val('1.1.2'), topic: versions, emit: versions_fastqstats


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastq-stats $args ${fastq} -x ${prefix}_per_base.tsv > ${prefix}_general.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}_general.tsv
    touch ${prefix}_per_base.tsv
    """
}
