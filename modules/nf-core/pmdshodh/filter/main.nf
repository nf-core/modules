// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process PMDSHODH_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8b/8bdc1b50afea6c4c60f579807586e2d1603ebf5bc06bc414ed116dfad9e8b8d3/data'
        : 'community.wave.seqera.io/library/pmdshodh:0.3.3--4d662c161631c5ad'}"

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
