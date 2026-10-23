// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FGUMI_EXTRACT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/74/74c00516b0152fbf7f721f7b8983024c652639d47c5e225c0ad0f6ca00773788/data'
:         'community.wave.seqera.io/library/fgumi:0.8.0--4fd382d441fde4a1' }"

    input:
    tuple val(meta), path(reads), val(library)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    fgumi extract \\
        --inputs ${reads.join(' ')} \\
        --output ${prefix}.bam \\
        ${args} \\
        --sample ${prefix} \\
        --library "${library}"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
