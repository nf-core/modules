process SKANI_SEARCH {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/skani:0.2.2--ha6fb395_2':
        'biocontainers/skani:0.2.2--ha6fb395_2' }"

    input:
    tuple val(meta) , path(query)
    tuple val(meta2), path(sketch_dir)

    output:
    tuple val(meta), path("${prefix}.tsv") , emit: search
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    skani \\
        search \\
            ${query} \\
            -d ${sketch_dir} \\
            -t ${task.cpus} \\
            -o ${prefix}.tsv \\
            ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version 2>&1 | sed 's/^.*skani //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version 2>&1 | sed 's/^.*skani //; s/ .*\$//')
    END_VERSIONS
    """
}
