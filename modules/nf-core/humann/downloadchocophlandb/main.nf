process HUMANN_DOWNLOADCHOCOPHLANDB {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann:3.8--pyh7cba7a3_0':
        'biocontainers/humann:3.8--pyh7cba7a3_0' }"

    input:
    val chocophlan_db_version

    output:
    path("chocophlan")  , emit: chocophlan_db
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    humann_databases \\
        --download chocophlan \\
        ${chocophlan_db_version} \\
        . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$(echo \$(humann --version 2>&1 | sed 's/^.*humann //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    mkdir chocophlan

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$(echo \$(humann --version 2>&1 | sed 's/^.*humann //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
