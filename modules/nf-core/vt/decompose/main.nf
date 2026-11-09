// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VT_DECOMPOSE {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c17ccba14c98db0059baa135debcea14a2669c18a67dbdef98fa7ba547f8344/data'
        : 'community.wave.seqera.io/library/htslib_vt:933758f4e7c3b9f1'}"

    input:
    tuple val(meta), path(vcf), path(intervals)

    output:
    tuple val(meta), path("*.vcf.gz")   , emit: vcf
    tuple val("${task.process}"), val('vt'), val('0.57721'), topic: versions, emit: versions_vt // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_decompose"

    if ("$vcf" == "${prefix}.vcf" || "$vcf" == "${prefix}.vcf.gz") {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    def bed = intervals ? "-i ${intervals}" : ""

    """
    vt decompose \\
        -o ${prefix}.vcf \\
        ${bed} \\
        ${args} \\
        ${vcf}

    bgzip ${args2} --threads ${task.cpus} ${prefix}.vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_decompose"

    if ("$vcf" == "${prefix}.vcf" || "$vcf" == "${prefix}.vcf.gz") {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    """
    echo "" | gzip > ${prefix}.vcf.gz
    """
}
