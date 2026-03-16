process SAMTOOLS_MPILEUP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e5/e5598451c6d348cce36191bafe1911ad71e440137d7a329da946f2b0dbb0e7f3/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23--cde2c40a51d6f752'}"

    input:
    tuple val(meta), path(input), path(index), path(intervals)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.mpileup.gz"), emit: mpileup
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_cmd = fasta ? "--fasta-ref ${fasta}" : ""
    def intervals_cmd = intervals ? "-l ${intervals}" : ""
    """
    samtools mpileup \\
        ${fasta_cmd} \\
        --output ${prefix}.mpileup \\
        ${args} \\
        ${intervals_cmd} \\
        ${input}
    bgzip ${prefix}.mpileup
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.mpileup.gz
    """
}
