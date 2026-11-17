// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

nextflow.enable.types = true

process PRESEQ_LCEXTRAP {
    tag "${id}"
    label 'process_single'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/49/497b9b64f7a0554e9e0eaba8b789fc8796a01f61547a52d6006774bea4b853d0/data'
        : 'community.wave.seqera.io/library/preseq:3.2.0--c485d1f68d51de10'}"

    ext [
        prefix: { output_prefix ?: "${id}" },
        args: { preseq_lcextrap_args ?: '' }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        single_end: Boolean,
        bam: Path,
        preseq_lcextrap_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        lc_extrap: file("*.lc_extrap.txt"),
        log: file("*.log")
    )

    topic:
    record(process: "${task.process}", tool: 'preseq', version: eval("preseq 2>&1 | sed -n 's/Version: //p'")) >> 'versions'

    script:
    // A retry disables testing for defects in the count histogram
    def args = task.attempt > 1 ? "${task.ext.args} -defects" : task.ext.args
    def paired_end = single_end ? '' : '-pe'
    """
    preseq \\
        lc_extrap \\
        ${args} \\
        ${paired_end} \\
        -output ${task.ext.prefix}.lc_extrap.txt \\
        ${bam}
    cp .command.err ${task.ext.prefix}.command.log
    """

    stub:
    """
    touch ${task.ext.prefix}.lc_extrap.txt
    touch ${task.ext.prefix}.command.log
    """
}
