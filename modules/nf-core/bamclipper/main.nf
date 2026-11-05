// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BAMCLIPPER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c5/c55c19c5134421e0c5c5dd846f30ca06074faef6c011806be7d35b94e19916bb/data'
        : 'community.wave.seqera.io/library/bamclipper_samtools:34439705a8946c0f' }"

    input:
    tuple val(meta), path(bam), path(bai), path(bedpe)

    output:
    tuple val(meta), path("*.primerclipped.bam")    , emit: bam
    tuple val(meta), path("*.primerclipped.bam.bai"), emit: bai
    tuple val("${task.process}"), val('bamclipper'), val("1.0.0"), emit: versions_bamclipper, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions
    // WARN: Version information not provided by tool on CLI
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bamclipper.sh \\
        -b ${bam} \\
        -p ${bedpe} \\
        -n $task.cpus \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.primerclipped.bam
    touch ${prefix}.primerclipped.bam.bai
    """
}
