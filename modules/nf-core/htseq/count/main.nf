process HTSEQ_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htseq:2.0.3--py310ha14a713_0':
        'biocontainers/htseq:2.0.3--py310ha14a713_0' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

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


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htseq: \$(echo \$(htseq-count --version ) | sed 's/^.*htseq-count //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(htseq-count --version ) | sed 's/^.*htseq-count //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
