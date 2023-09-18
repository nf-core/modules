
process UPP_ALIGN {
    tag "$meta_fasta.id"
    label 'process_medium'

    conda "bioconda::sepp bioconda::pasta"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta_fasta), path(fasta)
    tuple val(meta_tree), path(tree)

    output:
    tuple val(meta), path("*_alignment.fasta"), emit: alignment
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    meta = meta_tree + meta_fasta
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def loadtree = tree ? "-t $tree" : ''
    """
    python $(which sepp)/run_upp.py $args \\
        -x $task.cpus \\
        -o ${prefix}.aln \\
        $loadtree \\
        -s $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    touch ${prefix}.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
