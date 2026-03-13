process CONCOCT_CUTUPFASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/concoct:1.1.0--py39h8907335_8'
        : 'biocontainers/concoct:1.1.0--py39h8907335_8'}"

    input:
    tuple val(meta), path(fasta)
    val bed

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val(meta), path("*.bed"), optional: true, emit: bed
    tuple val("${task.process}"), val('concoct'), eval('concoct --version | cut -d " " -f2'), emit: versions_concoct, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bedfile = bed ? "-b ${prefix}.bed" : ""
    if ("${fasta}" == "${prefix}.fasta") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    cut_up_fasta.py \\
        ${fasta} \\
        ${args} \\
        ${bedfile} \\
        > ${prefix}.fasta

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${fasta}" == "${prefix}.fasta") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    touch ${prefix}.fasta

    """
}
