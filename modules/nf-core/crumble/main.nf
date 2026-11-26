// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CRUMBLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0a/0a1346c2b14945ca6fdd5f0ad16793359bda5515adb73ec9bc8fba29667f1c15/data'
        : 'community.wave.seqera.io/library/crumble_htslib:7794eee7eca9bbfb'}"

    input:
    tuple val(meta), path(input)
    path keepbed
    val bedout

    output:
    tuple val(meta), path("*.bam"),  emit: bam,     optional: true
    tuple val(meta), path("*.cram"), emit: cram,    optional: true
    tuple val(meta), path("*.sam"),  emit: sam,     optional: true
    tuple val(meta), path("*.bed"),  emit: bed,     optional: true
    //WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('crumble'), val('0.9.1'), emit: versions_crumble, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def extension   = args.contains("-O sam") ? "sam" :
                     args.contains("-O bam") ? "bam" :
                     args.contains("-O cram") ? "cram" :
                     "bam"
    def bedin       = keepbed ? "-R ${keepbed}" : ""
    def bedout_args = bedout ? "-b ${prefix}.out.bed" : ""
    if ("$input" == "${prefix}.${extension}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    crumble \\
        $args \\
        $bedin \\
        $bedout_args \\
        $input \\
        ${prefix}.${extension}
    """

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def extension   = args.contains("-O sam") ? "sam" :
                     args.contains("-O bam") ? "bam" :
                     args.contains("-O cram") ? "cram" :
                     "bam"
    def bedout_args = bedout ? "touch ${prefix}.out.bed" : ''
    if ("$input" == "${prefix}.${extension}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    touch ${prefix}.${extension}
    $bedout_args
    """
}
