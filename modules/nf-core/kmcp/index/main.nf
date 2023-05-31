process KMCP_INDEX {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kmcp=0.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kmcp:0.9.1--h9ee0642_0':
        'quay.io/biocontainers/kmcp:0.9.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(compute_dir)

    output:
    tuple val(meta), path("${prefix}")   , emit: kmcp
    tuple val(meta), path("*.log")       , emit: log
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    kmcp \\
        index \\
        --in-dir $compute_dir \\
        $args \\
        --threads $task.cpus \\
        --log ${prefix}.log \\
        --out-dir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmcp: \$(echo \$(kmcp version 2>&1) | sed -n 1p | sed 's/^.*kmcp v//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmcp: \$(echo \$(kmcp version 2>&1) | sed -n 1p | sed 's/^.*kmcp v//')
    END_VERSIONS
    """
}
