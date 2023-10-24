
process KMCP_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kmcp=0.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kmcp:0.9.1--h9ee0642_0':
        'quay.io/biocontainers/kmcp:0.9.1--h9ee0642_0' }"


    input:
    tuple val(meta), path(search_out)

    output:
    tuple val(meta), path("*.gz"),  emit: result
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kmcp \\
        merge \\
        $args \\
        --threads $task.cpus \\
        --out-file ${prefix}.gz \\
        $search_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmcp: \$(echo \$(kmcp version 2>&1) | sed -n 1p | sed 's/^.*kmcp v//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}
    gzip ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmcp: \$(echo \$(kmcp version 2>&1) | sed -n 1p | sed 's/^.*kmcp v//')
    END_VERSIONS
    """
}
