// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BEDTOOLS_BAMTOBED {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52d3fbb86c48b8c39a989733ca97bfd84cf665a9ab279be7cfbbb36f5a53c46a/data'
        : 'community.wave.seqera.io/library/bedtools:2.31.1--897920cd303d6413' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedtools \\
        bamtobed \\
        ${args} \\
        -i ${bam} \\
        > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    """
}
