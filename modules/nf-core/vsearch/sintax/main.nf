process VSEARCH_SINTAX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.21.1--h95f258a_0':
        'biocontainers/vsearch:2.21.1--h95f258a_0' }"

    input:
    tuple val(meta), path(queryfasta)
    path db

    output:
    tuple val(meta), path('*.tsv')   , optional: true, emit: tsv
    tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | sed -n "1s/.*v\\([0-9.]*\\).*/\\\\1/p"'), emit: versions_vsearch, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vsearch \\
        --sintax $queryfasta \\
        --db $db \\
        --threads $task.cpus \\
        $args \\
        --tabbedout ${prefix}.tsv
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
