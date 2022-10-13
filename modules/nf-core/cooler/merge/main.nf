process COOLER_MERGE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::cooler=0.8.11" : null)
    def container_image = "/cooler:0.8.11--pyh3252c3a_0"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(cool)

    output:
    tuple val(meta), path("*.cool"), emit: cool
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cooler merge \\
        $args \\
        ${prefix}.cool \\
        ${cool}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(cooler --version 2>&1 | sed 's/cooler, version //')
    END_VERSIONS
    """
}
