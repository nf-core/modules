def VERSION='2.7.1' // Version information not provided by tool on CLI

process KRONA_KTUPDATETAXONOMY {
    label 'process_single'

    conda (params.enable_conda ? "bioconda::krona=2.7.1" : null)
        def container_image = "/krona:2.7.1--pl526_5"
                                                         container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    output:
    path 'taxonomy/taxonomy.tab', emit: db
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ktUpdateTaxonomy.sh \\
        $args \\
        taxonomy/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krona: $VERSION
    END_VERSIONS
    """
}
