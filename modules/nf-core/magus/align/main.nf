process MAGUS_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ae4ea1182e75371808710b6c081bef8b228c4815:10b41722a6b9471a0945fe6baeb9aff444d8eb1d-0':
        'biocontainers/mulled-v2-ae4ea1182e75371808710b6c081bef8b228c4815:10b41722a6b9471a0945fe6baeb9aff444d8eb1d-0' }"

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
    def prefix = task.ext.prefix ?: "${meta.id}"
    def loadtree = tree ? "-t $tree" : ''
    def write_output = compress ? "--overwrite -o >(pigz -cp ${task.cpus} > ${prefix}.aln.gz)" : "-o ${prefix}.aln"
    // using >() is necessary to preserve the return value,
    // so nextflow knows to display an error when it failed
    // using --overwrite is necessary, as the file descriptor generated by the named file will alerady exist
    """
    magus \\
        -np $task.cpus \\
        -i $fasta \\
        -d ./ \\
        $write_output \\
        $loadtree \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MAGUS: \$(magus --version)
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln${compress ? '.gz' : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MAGUS: \$(magus --version)
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
