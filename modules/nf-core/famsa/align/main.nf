process FAMSA_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/famsa:2.4.1--h9ee0642_0':
        'biocontainers/famsa:2.4.1--h9ee0642_0' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(tree)
    val(compress)

    output:
    tuple val(meta), path("*.aln{.gz,}"), emit: alignment
    tuple val("${task.process}"), val('famsa'), eval("famsa -help 2>&1 | head -n 2 | tail -n 1 | sed 's/ version //' | awk '{print \$1}' | xargs"), topic: versions, emit: versions_famsa

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
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    ${compress ? "echo '' | gzip > ${prefix}.aln.gz" : "touch ${prefix}.aln"}
    """
}
