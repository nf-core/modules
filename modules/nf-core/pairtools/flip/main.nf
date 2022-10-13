process PAIRTOOLS_FLIP {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pairtools=0.3.0" : null)
        def container_image = "/pairtools:0.3.0--py37hb9c2fc3_5"
                                                 container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(sam)
    path chromsizes

    output:
    tuple val(meta), path("*.flip.gz"), emit: flip
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pairtools \\
        flip \\
        -c $chromsizes \\
        $args \\
        -o ${prefix}.flip.gz \\
        $sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version 2>&1 | sed 's/pairtools.*version //')
    END_VERSIONS
    """
}
