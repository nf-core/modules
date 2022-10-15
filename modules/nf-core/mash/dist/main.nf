process MASH_DIST {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mash=2.3" : null)
    def container_image = "mash:2.3--he348c14_1"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(query)
    path reference

    output:
    tuple val(meta), path("*.txt"), emit: dist
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mash \\
        dist \\
        -p $task.cpus \\
        $args \\
        $reference \\
        $query > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """
}
