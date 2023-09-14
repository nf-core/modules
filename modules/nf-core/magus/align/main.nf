process MAGUS_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::magus-msa=0.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magus-msa:0.1.0':
        'biocontainers/magus-msa:0.1.0' }"

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
    """
    #TODO make proper call, dry run on cluster
    magus \\
        -@ $task.cpus \\
        -i $ \\
        -o ${prefix}.aln \\
        -T $prefix \\
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
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
