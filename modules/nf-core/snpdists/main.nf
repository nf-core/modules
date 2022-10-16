process SNPDISTS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::snp-dists=0.8.2" : null)
    def container_image = "snp-dists:0.8.2--h5bf99c6_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    snp-dists \\
        $args \\
        $alignment > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpdists: \$(snp-dists -v 2>&1 | sed 's/snp-dists //;')
    END_VERSIONS
    """
}
