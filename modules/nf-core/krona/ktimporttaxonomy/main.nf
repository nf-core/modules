process KRONA_KTIMPORTTAXONOMY {
    tag "${meta.id}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::krona=2.8" : null)
    def container_image = "/krona:2.8--pl5262hdfd78af_2"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(report)
    path taxonomy, stageAs: 'taxonomy.tab'

    output:
    tuple val(meta), path ('*.html'), emit: html
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.8' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    TAXONOMY=\$(find -L . -name '*.tab' -exec dirname {} \\;)
    echo \$TAXONOMY

    ktImportTaxonomy \\
        $args \\
        -o ${prefix}.html \\
        -tax \$TAXONOMY/ \\
        $report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krona: $VERSION
    END_VERSIONS
    """
}
