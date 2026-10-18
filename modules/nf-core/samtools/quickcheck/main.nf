process SAMTOOLS_QUICKCHECK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/462b5db8b8a4a1a76a627a056a56a3c9ab344d51db27dee8fc89109d23bfea71/data'
        : 'community.wave.seqera.io/library/samtools_htslib:1.24--e186722e17fa9b47'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), env("EXIT_CODE"), emit: bam
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools quickcheck \\
        ${args} \\
        ${bam} || EXIT_CODE=\$?

    EXIT_CODE=\${EXIT_CODE:-0}
    """

    stub:
    def args = task.ext.args ?: ''
    """
    EXIT_CODE=0
    echo ${args}
    """
}
