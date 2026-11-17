// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process HTSNIMTOOLS_VCFCHECK {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2e/2e83030e0fd490b217e51050658481e37370e4b2e6dc68576467eb504626708d/data'
        : 'community.wave.seqera.io/library/hts-nim-tools:0.3.31--e705f8578bf971e8'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(background_vcf), path(background_tbi)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('htsnimtools'), eval("hts_nim_tools | sed -n 's/version: //p'"), emit: versions_htsnimtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hts_nim_tools \\
        vcf-check \\
        ${args} \\
        ${background_vcf} \\
        ${vcf} \\
        > ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
