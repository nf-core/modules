process CMSEQ_POLYMUT {
    tag "${meta.id}"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cmseq:1.0.4--pyhb7b1952_0'
        : 'biocontainers/cmseq:1.0.4--pyhb7b1952_0'}"

    input:
    tuple val(meta), path(bam), path(bai), path(gff), path(fasta)

    output:
    tuple val(meta), path("*.txt"), emit: polymut
    tuple val("${task.process}"), val('cmseq'), val(VERSION), topic: versions, emit: versions_cmseq

    when:
    task.ext.when == null || task.ext.when

    script:
    VERSION = '1.0.4'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_refid = fasta ? "-c ${fasta}" : ""
    def sortindex = bai ? "" : "--sortindex"
    """
    polymut.py \\
        ${args} \\
        ${sortindex} \\
        ${fasta_refid} \\
        --gff_file ${gff} \\
        ${bam} > ${prefix}.txt
    """

    stub:
    VERSION = '1.0.4'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
