// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ADMIXTURE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/71/714f67e96b6476475df13d17d2f1b69c5e7fecf04bee6da07e45776e165587b7/data'
        : 'community.wave.seqera.io/library/admixture:1.3.0--a0588fc0dcd0bd38'}"

    input:
    tuple val(meta), path (bed_ped_geno), path(bim_map), path(fam)
    val K


    output:
    tuple val(meta), path("*.Q"), emit: ancestry_fractions
    tuple val(meta), path("*.P"), emit: allele_frequencies
    tuple val("${task.process}"), val('admixture'), eval('admixture --version | tail -n 1'), emit: versions_admixture, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    """
    admixture \\
        ${bed_ped_geno} \\
        ${K} \\
        -j${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.Q"
    touch "${prefix}.P"
    """
}
