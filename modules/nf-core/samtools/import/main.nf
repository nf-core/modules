// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SAMTOOLS_IMPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/87660348a4d032cc08d6ac8101ff9d4d27164300d39c9472043bb5001f9ba801/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.25--704fb80d75745bd4'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.sam"), emit: sam, optional: true
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("--output-fmt sam")
        ? "sam"
        : args.contains("--output-fmt bam")
            ? "bam"
            : args.contains("--output-fmt cram")
                ? "cram"
                : "bam"
    def input = reads instanceof List && meta.single_end
        ? reads.join(" -0")
        : reads instanceof List && !meta.single_end
            ? "-1 ${reads[0]} -2 ${reads[1]}"
            : meta.single_end
                ? "-0 ${reads}"
                : !meta.single_end
                    ? "-s ${reads}"
                    : reads
    // if all else fails, just add the reads without flags
    """
    samtools \\
        import \\
        ${input} \\
        ${args} \\
        -@ ${task.cpus} \\
        -o ${prefix}.${suffix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bam
    """
}
