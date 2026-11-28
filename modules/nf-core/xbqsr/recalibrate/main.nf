// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process XBQSR_RECALIBRATE {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'

    container "ghcr.io/ikabadzhov/xbqsr:cuda"

    input:
    tuple val(meta), path(input), path(input_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(known_sites)

    output:
    tuple val(meta), path("*.table"), emit: table
    tuple val("${task.process}"), val('xbqsr'), eval("bqsr --version | sed 's/bqsr //'"), topic: versions, emit: versions_xbqsr

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("X-BQSR module does not support Conda. Please use Docker / Singularity / Apptainer with -profile gpu instead.")
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sites  = known_sites instanceof List ? known_sites.join(',') : "${known_sites}"
    """
    export HOME=\$PWD
    export XDG_DATA_HOME=\$PWD/.cache

    bqsr recalibrate \\
        --bam-in ${input} \\
        --recal-out ${prefix}.table \\
        --ref ${fasta} \\
        --known-sites ${sites} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.table
    """
}
