process HTSEQ_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htseq:2.0.3--py310ha14a713_0':
        'quay.io/biocontainers/htseq:2.0.3--py310ha14a713_0' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('htseq'), eval("htseq-count --version"), emit: versions_htseq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    htseq-count \\
        ${input} \\
        ${gtf} \\
        ${args} \\
        > ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
