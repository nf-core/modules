// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process XBQSR_APPLY {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'

    container "ghcr.io/ikabadzhov/xbqsr:cuda"

    input:
    tuple val(meta), path(input), path(input_index), path(table)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('xbqsr'), eval("bqsr --version | sed 's/bqsr //'"), topic: versions, emit: versions_xbqsr

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("X-BQSR module does not support Conda. Please use Docker / Singularity / Apptainer with -profile gpu instead.")
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.recal"
    if ("${input}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    export HOME=\$PWD
    export XDG_DATA_HOME=\$PWD/.cache

    bqsr apply \\
        --bam-in ${input} \\
        --bam-out ${prefix}.bam \\
        --recal ${table} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.recal"
    """
    touch ${prefix}.bam
    """
}
