// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ALDY_GENOTYPE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e8/e8ad9095fa0372544cb2c9d868620c5c5226b7c6504fa3c59a879d73bef52e1d/data'
        : 'community.wave.seqera.io/library/mappy_pysam_gcc_pip_pruned:76bc0744bd24414b'}"

    input:
    tuple val(meta), path(input), path(index), val(gene), val(profile)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.aldy"), emit: aldy
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('aldy'), eval("aldy help 2>&1 | grep -o 'Aldy v[0-9.]*' | cut -c7-"), topic: versions, emit: versions_aldy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ''
    """
    aldy \\
        genotype \\
        ${args} \\
        --gene ${gene} \\
        --profile ${profile} \\
        ${reference} \\
        --output ${prefix}.aldy \\
        --log ${prefix}.log \\
        ${input}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.aldy
    touch ${prefix}.log
    """
}
