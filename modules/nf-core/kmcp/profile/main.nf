process KMCP_PROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kmcp=0.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kmcp:0.9.1--h9ee0642_0':
        'biocontainers/kmcp:0.9.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(search_results)
    path (db)
    val mode

    output:
    tuple val(meta), path("*.profile"), emit: profile
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    taxid=`find -L ${db} -name "*map"`
    taxdump=`find -L ${db}/*/ -type d  -not -name "R001"`
    kmcp \\
        profile \\
        $args \\
        -X \$taxdump \\
        -T \$taxid \\
        -m $mode \\
        -j $task.cpus \\
        -o ${prefix}.profile \\
        $search_results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmcp: \$(echo \$(kmcp version 2>&1) | sed -n 1p | sed 's/^.*kmcp v//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.profile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmcp: \$(echo \$(kmcp version 2>&1) | sed -n 1p | sed 's/^.*kmcp v//')
    END_VERSIONS
    """
}
