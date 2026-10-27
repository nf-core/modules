// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VT_DECOMPOSEBLOCKSUB {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/15/15c33ba0e0d0e953d83fa810e10ebb0e1617a539ff82660ba7b641154fd3d25b/data'
        : 'community.wave.seqera.io/library/vt:0.57721--83c629c7e3923326'}"

    input:
    tuple val(meta), path(vcf), path(index), path(intervals)

    output:
    tuple val(meta), path("*.vcf.gz")   , emit: vcf
    tuple val("${task.process}"), val('vt'), val('0.57721'), topic: versions, emit: versions_vt // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_decomposeblocksub"

    if ("$vcf" == "${prefix}.vcf.gz") {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    def bed = intervals ? "-i ${intervals}" : ""
    """
    vt decompose_blocksub \\
        -o ${prefix}.vcf.gz \\
        ${bed} \\
        ${args} \\
        ${vcf}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_decomposeblocksub"

    if ("$vcf" == "${prefix}.vcf.gz") {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    """
    echo "" | gzip > ${prefix}.vcf.gz
    """
}
