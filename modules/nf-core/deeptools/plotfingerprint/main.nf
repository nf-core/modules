// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process DEEPTOOLS_PLOTFINGERPRINT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ba/ba9641d61acf134c6b989b3af5f25bcd5aadb3faf692a41c00ba4cd25f99032e/data'
:         'community.wave.seqera.io/library/deeptools:3.5.6--511a2b05870b9ec3' }"

    input:
    tuple val(meta), path(bams), path(bais)

    output:
    tuple val(meta), path("*.pdf")          , emit: pdf
    tuple val(meta), path("*.raw.txt")      , emit: matrix
    tuple val(meta), path("*.qcmetrics.txt"), emit: metrics
    tuple val("${task.process}"), val('deeptools'), eval('plotFingerprint --version | sed "s/plotFingerprint //g"') , emit: versions_deeptools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extend   = (meta.single_end && params.fragment_size > 0) ? "--extendReads ${params.fragment_size}" : ''
    """
    plotFingerprint \\
        $args \\
        $extend \\
        --bamfiles ${bams.join(' ')} \\
        --plotFile ${prefix}.plotFingerprint.pdf \\
        --outRawCounts ${prefix}.plotFingerprint.raw.txt \\
        --outQualityMetrics ${prefix}.plotFingerprint.qcmetrics.txt \\
        --numberOfProcessors $task.cpus
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.plotFingerprint.pdf
    touch ${prefix}.plotFingerprint.raw.txt
    touch ${prefix}.plotFingerprint.qcmetrics.txt
    """
}
