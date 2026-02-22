process SEQFU_CHECK {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqfu:1.22.3--hc29b5fc_1'
        : 'biocontainers/seqfu:1.22.3--hc29b5fc_1'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: check
    tuple val("${task.process}"), val('seqfu'), eval('seqfu version'), emit: versions_seqfu, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def dirFlag = (reads instanceof List ? reads.every { read -> read.isDirectory() } : reads.isDirectory()) ? "--dir" : ""
    """
    seqfu \\
        check \\
        ${args} \\
        ${dirFlag} ${reads} > ${prefix}.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
