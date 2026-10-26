// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f2/f26f1ce0aac66b8f199594c1ed6bb39f256b5a08d524167e84176fb5799f776a/data'
:         'community.wave.seqera.io/library/cutadapt:5.2--efafca34ac3f6c8a' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')          , emit: log
    tuple val("${task.process}"), val("cutadapt"), eval('cutadapt --version'), topic: versions, emit: versions_cutadapt

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed  = meta.single_end ? "-o ${prefix}.trim.fastq.gz" : "-o ${prefix}_1.trim.fastq.gz -p ${prefix}_2.trim.fastq.gz"
    """
    cutadapt \\
        --cores $task.cpus \\
        $args \\
        $trimmed \\
        $reads \\
        > ${prefix}.cutadapt.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        output_command = "echo '' | gzip > ${prefix}.trim.fastq.gz ;"
    }
    else {
        output_command  = "echo '' | gzip > ${prefix}_1.trim.fastq.gz ;"
        output_command += "echo '' | gzip > ${prefix}_2.trim.fastq.gz ;"
    }
    """
    ${output_command}
    touch ${prefix}.cutadapt.log
    """
}
