process SKANI_TRIANGLE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/skani:0.2.1--h4ac6f70_0':
        'biocontainers/skani:0.2.1--h4ac6f70_0' }"

    input:
    tuple val(meta) , path(queries)

    output:
    tuple val(meta), path("${prefix}.tsv")  , emit: triangle
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    skani \\
        triangle \\
            ${queries} \\
            --sparse \\
            -o ${prefix}.tsv \\
            -t ${task.cpus} \\
            ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version 2>&1 | sed 's/^.*skani //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version 2>&1 | sed 's/^.*skani //; s/ .*\$//')
    END_VERSIONS
    """
}
