process KMCP_SEARCH {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kmcp=0.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kmcp:0.9.1--h9ee0642_0':
        'biocontainers/kmcp:0.9.1--h9ee0642_0' }"

    input:
    path(db)
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.gz") , emit: result
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input  = meta.single_end ? "${reads}": "-1 ${reads[0]} -2 ${reads[1]}"
    """
    kmcp \\
        search \\
        $args \\
        --threads $task.cpus \\
        --db-dir $db \\
        $reads \\
        --out-file ${prefix}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmcp: \$(echo \$(kmcp version 2>&1) | sed -n 1p | sed 's/^.*kmcp v//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmcp: \$(echo \$(kmcp version 2>&1) | sed -n 1p | sed 's/^.*kmcp v//')
    END_VERSIONS
    """

}
