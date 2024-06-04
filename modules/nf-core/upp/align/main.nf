process UPP_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sepp:4.5.1--py38h9ee0642_4':
        'biocontainers/sepp:4.5.1--py38h9ee0642_4' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(tree)

    output:
    tuple val(meta), path("*.aln"), emit: alignment
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
