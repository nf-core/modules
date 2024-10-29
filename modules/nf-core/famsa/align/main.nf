

process FAMSA_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/famsa:2.2.2--h9f5acd7_0':
        'biocontainers/famsa:2.2.2--h9f5acd7_0' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(tree)
    val(compress)

    output:
    tuple val(meta), path("*.aln{.gz,}"), emit: alignment
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def compress_args = compress ? '-gz' : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def options_tree = tree ? "-gt import $tree" : ""
    """
    famsa $options_tree \\
        $compress_args \\
        $args \\
        -t ${task.cpus} \\
        ${fasta} \\
        ${prefix}.aln${compress ? '.gz':''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        famsa: \$( famsa -help 2>&1 | head -n 2 | tail -n 1 | sed 's/ version //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln${compress ? '.gz' : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        famsa: \$( famsa -help 2>&1 | head -n 2 | tail -n 1 | sed 's/ version //g' )
    END_VERSIONS
    """
}
