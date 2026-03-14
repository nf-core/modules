process SAMTOOLS_BAM2FQ {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e5/e5598451c6d348cce36191bafe1911ad71e440137d7a329da946f2b0dbb0e7f3/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23--cde2c40a51d6f752'}"

    input:
    tuple val(meta), path(inputbam)
    val split

    output:
    tuple val(meta), path("*.fq.gz"), emit: reads
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

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
            ${inputbam} | gzip --no-name > ${prefix}_interleaved.fq.gz
        """
    }

    stub:

    def create_cmd = "echo | gzip >"
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (split) {
        """
        ${create_cmd} ${prefix}_1.fq.gz
        ${create_cmd} ${prefix}_2.fq.gz
        ${create_cmd} ${prefix}_other.fq.gz
        ${create_cmd} ${prefix}_singleton.fq.gz
        """
    }
    else {
        """
        ${create_cmd} ${prefix}_interleaved.fq.gz
        """
    }
}
