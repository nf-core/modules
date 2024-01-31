process QUARTONOTEBOOK {
    tag "$meta.id"
    label 'process_low'

    // TODO: Harmonise Quarto versions between Docker/Conda
    conda "conda-forge::quarto=1.3.433 conda-forge::r-base=4.3.2 conda-forge::r-rmarkdown=2.25 conda-forge::matplotlib=3.4.3"
    container "docker.io/erikfas/quartonotebook"

    input:
    tuple val(meta), path(notebook)
    path input_files

    output:
    tuple val(meta), path("*.html")     , emit: html
    tuple val(meta), path("artifacts/*"), emit: artifacts, optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    quarto render ${notebook}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto -v)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto -v)
    END_VERSIONS
    """
}
