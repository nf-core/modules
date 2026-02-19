
process VSEARCH_FASTQFILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.28.1--h6a68c12_1':
        'biocontainers/vsearch:2.28.1--h6a68c12_1' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path('*.fasta')   , emit: fasta
    path "*.log"                       , emit: log
    tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | sed -n "1s/.*v\\([0-9.]*\\).*/\\\\1/p"'), emit: versions_vsearch, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vsearch \\
        --fastq_filter ${fastq} \\
        $args \\
        --fastaout ${prefix}.fasta 2>&1 | tee ${prefix}.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta
    touch ${prefix}.log
    """
}
