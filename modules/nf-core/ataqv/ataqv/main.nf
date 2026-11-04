// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ATAQV_ATAQV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f7063be8c7419222c75482548d138078837657f628d149e1a0b3e6a9080799de/data'
        : 'community.wave.seqera.io/library/ataqv:1.3.1--13b9c1b1eb8ada3f'}"

    input:
    tuple val(meta), path(bam), path(bai), path(peak_file)
    val organism
    val mito_name
    path tss_file
    path excl_regs_file
    path autosom_ref_file

    output:
    tuple val(meta), path("*.ataqv.json"), emit: json
    tuple val(meta), path("*.problems")  , emit: problems, optional: true
    tuple val("${task.process}"), val('ataqv'), eval("ataqv --version 2>&1 || true"), emit: versions_ataqv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mito   = mito_name ? "--mitochondrial-reference-name ${mito_name}" : ''

    def peak        = peak_file        ? "--peak-file ${peak_file}"                       : ''
    def tss         = tss_file         ? "--tss-file ${tss_file}"                         : ''
    def excl_regs   = excl_regs_file   ? "--excluded-region-file ${excl_regs_file}"       : ''
    def autosom_ref = autosom_ref_file ? "--autosomal-reference-file ${autosom_ref_file}" : ''
    """
    ataqv \\
        ${args} \\
        ${mito} \\
        ${peak} \\
        ${tss} \\
        ${excl_regs} \\
        ${autosom_ref} \\
        --metrics-file "${prefix}.ataqv.json" \\
        --threads ${task.cpus} \\
        --name ${prefix} \\
        ${organism} \\
        ${bam}

    """

    stub:
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def args          = task.ext.args   ?: ''
    def problems_cmd  = args.contains("--log-problematic-reads") ? "touch 1.problems" : ""
    """
    touch ${prefix}.ataqv.json
    ${problems_cmd}

    """
}
