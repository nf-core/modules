process LEGSTA {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::legsta=0.5.1" : null)
    def container_image = "/legsta:0.5.1--hdfd78af_2"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    legsta \\
        $args \\
        $seqs > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        legsta: \$(echo \$(legsta --version 2>&1) | sed 's/^.*legsta //; s/ .*\$//;')
    END_VERSIONS
    """
}
