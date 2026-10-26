// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SAMTOOLS_BAM2FQ {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/39db574ae6434e35bfb06db773f8bb7496a8109e7cc49bedc9e9947b3c87c797/data'
        : 'community.wave.seqera.io/library/htslib_samtools_xz:63148285647cb8e3'}"

    input:
    tuple val(meta), path(inputbam)
    val split

    output:
    tuple val(meta), path("*.fq.gz"), emit: reads
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools
    tuple val("${task.process}"), val("bgzip"), eval('bgzip --version | head -1 | sed "s/bgzip (htslib) //"'), emit: versions_bgzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (split) {
        """
        samtools \\
            bam2fq \\
            ${args} \\
            -@ ${task.cpus} \\
            -1 ${prefix}_1.fq.gz \\
            -2 ${prefix}_2.fq.gz \\
            -0 ${prefix}_other.fq.gz \\
            -s ${prefix}_singleton.fq.gz \\
            ${inputbam}
        """
    }
    else {
        """
        samtools \\
            bam2fq \\
            ${args} \\
            -@ ${task.cpus} \\
            ${inputbam} | bgzip > ${prefix}_interleaved.fq.gz
        """
    }

    stub:

    def prefix = task.ext.prefix ?: "${meta.id}"
    def bgzip_command_1 = "echo | bgzip -c > ${prefix}_1.fq.gz"
    def bgzip_command_2 = "echo | bgzip -c > ${prefix}_2.fq.gz"
    def bgzip_command_other = "echo | bgzip -c > ${prefix}_other.fq.gz"
    def bgzip_command_singleton = "echo | bgzip -c > ${prefix}_singleton.fq.gz"
    def bgzip_command_interleaved = "echo | bgzip -c > ${prefix}_interleaved.fq.gz"

    if (split) {
        """
        ${bgzip_command_1}
        ${bgzip_command_2}
        ${bgzip_command_other}
        ${bgzip_command_singleton}
        """
    }
    else {
        """
        ${bgzip_command_interleaved}
        """
    }
}
