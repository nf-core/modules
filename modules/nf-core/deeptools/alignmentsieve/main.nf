// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process DEEPTOOLS_ALIGNMENTSIEVE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/aa/aa2a707745965f561e777fdce54f18881d0e53a12c1b3a29ba1b2a6f2112f099/data'
:         'community.wave.seqera.io/library/deeptools_samtools:22ef5f45259c4336' }"

    input:
    tuple val(meta), path(input), path(input_index)

    output:
    tuple val(meta), path("*_as.bam") , emit: bam
    tuple val("${task.process}"), val('deeptools'), eval('alignmentSieve --version | sed "s/alignmentSieve //g"') , emit: versions_deeptools, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'") , emit: versions_samtools, topic: versions
    path  "*_log.txt"                 , emit: logs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    alignmentSieve \\
        $args \\
        -b $input \\
        -o ${prefix}_as.bam \\
        --filterMetrics ${prefix}_log.txt \\
        --numberOfProcessors $task.cpus
    """

    stub:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_as.bam
    touch ${prefix}_log.txt
    """
}
