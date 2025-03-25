process KMCP_COMPUTE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kmcp:0.9.4--h9ee0642_0':
        'biocontainers/kmcp:0.9.4--h9ee0642_0' }"

    input:
    tuple val(meta), path(sequences)

    output:
    tuple val(meta), path("${prefix}")             , emit: outdir
    tuple val(meta), path("${prefix}/_info.txt")   , emit: info
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input = sequences.isDirectory()? "--in-dir ${sequences}" : "${sequences}"
    """
    kmcp \\
        compute \\
        $args \\
        --threads $task.cpus \\
        --out-dir ${prefix}/ \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmcp: \$(echo \$(kmcp version 2>&1) | sed -n 1p | sed 's/^.*kmcp v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir $prefix
    touch $prefix/_info.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmcp: \$(echo \$(kmcp version 2>&1) | sed -n 1p | sed 's/^.*kmcp v//')
    END_VERSIONS
    """
}
