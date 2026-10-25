process SAMTOOLS_MPILEUP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/39db574ae6434e35bfb06db773f8bb7496a8109e7cc49bedc9e9947b3c87c797/data'
        : 'community.wave.seqera.io/library/htslib_samtools_xz:63148285647cb8e3'}"

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
    def bgzip_command = "echo | bgzip -c > ${prefix}.mpileup.gz"
    """
    ${bgzip_command}
    """
}
