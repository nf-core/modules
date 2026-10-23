// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process YAHS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/47e6425421c6823260bf4bd8d2d431f1f509e9bc9110bb3aa6b4652b85e63eac/data'
:         'community.wave.seqera.io/library/yahs:1.2.2--5a61538ab6dbd63d' }"

    input:
    tuple val(meta), path(fasta), path(fai), path(hic_map), path(agp)

    output:
    // note: typo in yahs file outputs - it writes "inital", not "initial"
    tuple val(meta), path("${prefix}_scaffolds_final.fa")    , emit: scaffolds_fasta   ,  optional: true
    tuple val(meta), path("${prefix}_scaffolds_final.agp")   , emit: scaffolds_agp     ,  optional: true
    tuple val(meta), path("${prefix}_{inital,no}_break*.agp"), emit: initial_break_agp ,  optional: true
    tuple val(meta), path("${prefix}_r*_*.agp")              , emit: round_agp         ,  optional: true
    tuple val(meta), path("${prefix}.bin")                   , emit: binary
    tuple val(meta), path("${prefix}.log")                   , emit: log
    tuple val("${task.process}"), val('yahs'), eval("yahs --version 2>&1"), emit: versions_yahs, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    def agp_input = agp ? "-a ${agp}" : ""
    """
    yahs \\
        -o ${prefix} \\
        ${agp_input} \\
        ${args} \\
        ${fasta} \\
        ${hic_map} \\
        2>| >( tee ${prefix}.log >&2 )
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_scaffolds_final.fa
    touch ${prefix}_scaffolds_final.agp
    touch ${prefix}_inital_break_01.agp
    touch ${prefix}_no_break.agp
    touch ${prefix}_r01.agp
    touch ${prefix}_r01_break.agp
    touch ${prefix}.bin
    touch ${prefix}.log
    """
}
