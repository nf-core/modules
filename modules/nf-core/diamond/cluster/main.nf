process DIAMOND_CLUSTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:2.1.9--h43eeafb_0':
        'biocontainers/diamond:2.1.9--h43eeafb_0' }"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem    = task.memory.toKilo() + 'K'
    def memarg = "-M ${mem}"
    """
    diamond \\
        cluster \\
        $args \\
        $memarg \\
        -p $task.cpus \\
        -d $db \\
        -o ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version |& sed '1!d ; s/diamond version //')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version |& sed '1!d ; s/diamond version //')
    END_VERSIONS
    """
}
