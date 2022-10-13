process MASHTREE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::mashtree=1.2.0" : null)
    def container_image = "/mashtree:1.2.0--pl526h516909a_0"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }
n

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("*.dnd"), emit: tree
    tuple val(meta), path("*.tsv"), emit: matrix
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mashtree \\
        $args \\
        --numcpus $task.cpus \\
        --outmatrix ${prefix}.tsv \\
        --outtree ${prefix}.dnd \\
        $seqs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mashtree: \$( echo \$( mashtree --version 2>&1 ) | sed 's/^.*Mashtree //' )
    END_VERSIONS
    """
}
