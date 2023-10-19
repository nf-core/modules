process UPP_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::sepp"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(tree)

    output:
    tuple val(meta), path("*_alignment.fasta"), emit: alignment
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def loadtree = tree ? "-t $tree" : ''
    """
    run_upp.py $args \\
        -x $task.cpus \\
        -o ${prefix}.aln \\
        $loadtree \\
        -s $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UPP: \$(echo \$(run_upp.py --version))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UPP: \$(echo \$(run_upp.py --version))
    END_VERSIONS
    """
}
