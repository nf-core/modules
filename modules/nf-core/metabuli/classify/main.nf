process METABULI_CLASSIFY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabuli:1.0.5--pl5321h6a68c12_1':
        'biocontainers/metabuli:1.0.5--pl5321h6a68c12_1' }"

    input:
    tuple val(meta), path(fastas)
    tuple val(meta2), path(db)

    output:
    tuple val(meta), path("*classifications.tsv"), emit: classification
    tuple val(meta), path("*report.tsv")         , emit: report
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    metabuli \\
        classify \\
        $args \\
        ${fastas} \\
        --threads $task.cpus \\
        --max-ram ${task.memory.toGiga()} \\
        ${db} \\
        . \\
        ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabuli: \$(metabuli version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_classifications.tsv"
    touch "${prefix}_report.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabuli: \$(metabuli version)
    END_VERSIONS
    """
}
