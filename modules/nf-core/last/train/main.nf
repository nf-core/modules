process LAST_TRAIN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/last:1542--h43eeafb_1' :
        'biocontainers/last:1542--h43eeafb_1' }"

    input:
    tuple val(meta), path(fastx)
    path  index

    output:
    tuple val(meta), path("*.train"), emit: param_file
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX_NAME=\$(basename \$(ls $index/*.des) .des)

    last-train \\
        $args \\
        -P $task.cpus \\
        ${index}/\$INDEX_NAME \\
        $fastx \\
        > ${prefix}.\$INDEX_NAME.train

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastdb --version | sed 's/lastdb //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX_NAME=STUB
    touch ${prefix}.\$INDEX_NAME.train

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastdb --version | sed 's/lastdb //')
    END_VERSIONS
    """
}
