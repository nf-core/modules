process MAGUS_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::magus-msa=0.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magus-msa:0.1.1':
        'biocontainers/magus-msa:0.1.1' }"

    input:
    tuple val(meta_fasta), path(fasta)
    tuple val(meta_tree), path(tree)

    output:
    tuple val(meta), path("*.aln"), emit: alignment
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def loadtree = tree ? "-t $tree" : ''
    """
    #TODO make proper call, dry run on cluster
    magus \\
        -np $task.cpus \\
        -i $fasta \\
        -o ${prefix}.aln \\
        $loadtree \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(magus --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(magus --version)
    END_VERSIONS
    """
}
