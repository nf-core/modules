// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process PRESEQ_LCEXTRAP {
    tag "$meta.id"
    label 'process_single'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/49/497b9b64f7a0554e9e0eaba8b789fc8796a01f61547a52d6006774bea4b853d0/data'
        : 'community.wave.seqera.io/library/preseq:3.2.0--c485d1f68d51de10'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.lc_extrap.txt"), emit: lc_extrap
    tuple val(meta), path("*.log")          , emit: log
    tuple val("${task.process}"), val('preseq'), eval("preseq 2>&1 | sed -n 's/Version: //p'"), emit: versions_preseq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    args = task.attempt > 1 ? args.join(' -defects') : args  // Disable testing for defects
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired_end = meta.single_end ? '' : '-pe'
    """
    preseq \\
        lc_extrap \\
        ${args} \\
        ${paired_end} \\
        -output ${prefix}.lc_extrap.txt \\
        ${bam}
    cp .command.err ${prefix}.command.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.lc_extrap.txt
    touch ${prefix}.command.log
    """
}
