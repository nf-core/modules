// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2a/2aac61eceef146ec8b1c0bb117542dfffb0dce0b720efb867c69a71ee4774f57/data'
        : 'community.wave.seqera.io/library/trimmomatic:0.39--1ffbebb3556a75f8'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.paired.trim*.fastq.gz")   , emit: trimmed_reads
    tuple val(meta), path("*.unpaired.trim_*.fastq.gz"), emit: unpaired_reads, optional:true
    tuple val(meta), path("*_trim.log")                , emit: trim_log
    tuple val(meta), path("*_out.log")                 , emit: out_log
    tuple val(meta), path("*.summary")                 , emit: summary
    tuple val("${task.process}"), val('trimmomatic'), eval("trimmomatic -version"), topic: versions, emit: versions_trimmomatic

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed = meta.single_end ? "SE" : "PE"
    def output = meta.single_end ?
        "${prefix}.SE.paired.trim.fastq.gz" // HACK to avoid unpaired and paired in the trimmed_reads output
        : "${prefix}.paired.trim_1.fastq.gz ${prefix}.unpaired.trim_1.fastq.gz ${prefix}.paired.trim_2.fastq.gz ${prefix}.unpaired.trim_2.fastq.gz"
    def qual_trim = task.ext.args2 ?: ''
    """
    trimmomatic \\
        $trimmed \\
        -threads $task.cpus \\
        -trimlog ${prefix}_trim.log \\
        -summary ${prefix}.summary \\
        $reads \\
        $output \\
        $qual_trim \\
        $args 2>| >(tee ${prefix}_out.log >&2)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        output_command = "echo '' | gzip > ${prefix}.SE.paired.trim.fastq.gz"
    } else {
        output_command  = "echo '' | gzip > ${prefix}.paired.trim_1.fastq.gz\n"
        output_command += "echo '' | gzip > ${prefix}.paired.trim_2.fastq.gz\n"
        output_command += "echo '' | gzip > ${prefix}.unpaired.trim_1.fastq.gz\n"
        output_command += "echo '' | gzip > ${prefix}.unpaired.trim_2.fastq.gz"
    }

    """
    $output_command
    touch ${prefix}.summary
    touch ${prefix}_trim.log
    touch ${prefix}_out.log
    """

}
