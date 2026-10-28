// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SAMTOOLS_CONVERT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c2f02ceb8cbe5aa37930dcbc0de6f1d07dfb292178c7e780c44ab2f556bf829f/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--8a3ddf262165ba23'}"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    tuple val(meta), path("*.crai"), emit: crai, optional: true
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_extension = input.getExtension() == "bam" ? "cram" : "bam"

    """
    samtools view \\
        --threads ${task.cpus} \\
        --reference ${fasta} \\
        ${args} \\
        ${input} \\
        -o ${prefix}.${output_extension}

    samtools index -@${task.cpus} ${prefix}.${output_extension}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_extension = input.getExtension() == "bam" ? "cram" : "bam"
    def index_extension = output_extension == "bam" ? "bai" : "crai"

    """
    touch ${prefix}.${output_extension}
    touch ${prefix}.${output_extension}.${index_extension}
    """
}
