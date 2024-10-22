process SOURCEPREDICT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourcepredict:0.5.1--pyhdfd78af_0':
        'biocontainers/sourcepredict:0.5.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(kraken_parse)
    path sources
    path labels

    output:
    tuple val(meta), path("*.sourcepredict.csv"), emit: report
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sourcepredict \\
        -s $sources \\
        -l $labels \\
        $args \\
        -t $task.cpus \\
        -o ${prefix}.sourcepredict.csv \\
        ${kraken_parse}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourcepredict: \$(echo \$(sourcepredict --help 2>&1 | head -n 1 | sed -n "s/.*SourcePredict v\\([0-9.]*\\).*/\\1/p")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.sourcepredict.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourcepredict: \$(echo \$(sourcepredict --help 2>&1 | head -n 1 | sed -n "s/.*SourcePredict v\\([0-9.]*\\).*/\\1/p")
    END_VERSIONS
    """
}
