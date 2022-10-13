process SEQKIT_STATS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::seqkit=2.2.0" : null)
        def container_image = "/seqkit:2.2.0--h9ee0642_0"
                                               container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.tsv"), emit: stats
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--all'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqkit stats \\
        --tabular \\
        $args \\
        $reads > '${prefix}.tsv'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}
