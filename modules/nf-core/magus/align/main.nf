process MAGUS_ALIGN {
    tag "$meta_fasta.id"
    label 'process_medium'

    conda "bioconda::magus-msa=0.1.2a"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magus-msa:0.1.2a':
        'biocontainers/magus-msa:0.1.2a' }"

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
    def prefix = task.ext.prefix ? "alignment" : "${meta_fasta.id}"
    def loadtree = tree ? "-t $tree" : ''
    meta = meta_tree + meta_fasta
    """
    magus \\
        -np $task.cpus \\
        -i $fasta \\
        -o ${prefix}.aln \\
        $loadtree \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MAGUS:\$(magus --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    meta = meta_tree + meta_fasta
    """
    touch ${prefix}.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MAGUS:\$(magus --version)
    END_VERSIONS
    """
}
