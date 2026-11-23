// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process PMDSHODH_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pmdshodh:0.3.1--hfa8f182_0'
        : 'quay.io/biocontainers/pmdshodh:0.3.1--hfa8f182_0'}"

    input:
    tuple val(meta), path(input)
    val threshold
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.json"), emit: json
    tuple val("${task.process}"), val('pmdshodh'), eval("pmdshodh --version | sed 's/^pmdshodh //'"), topic: versions, emit: versions_pmdshodh

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ''
    if ("${input}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    pmdshodh \\
        filter \\
        --threshold ${threshold} \\
        ${reference} \\
        --summary ${prefix}.pmdshodh.json \\
        --threads ${task.cpus} \\
        ${args} \\
        --output ${prefix}.bam \\
        ${input}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${input}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.bam
    touch ${prefix}.pmdshodh.json
    """
}
