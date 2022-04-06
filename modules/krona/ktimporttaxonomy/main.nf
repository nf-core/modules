def VERSION = '2.8' // Version information not provided by tool on CLI

process KRONA_KTIMPORTTAXONOMY {
    tag "${meta.id}"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::krona=2.8" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krona:2.8--pl5262hdfd78af_2' :
        'quay.io/biocontainers/krona:2.8--pl5262hdfd78af_2' }"

    input:
    tuple val(meta), path(report)
    path  "taxonomy/taxonomy.tab"

    output:
    tuple val(meta), path ('*.html'), emit: html
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ktImportTaxonomy "$report" -tax taxonomy

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krona: $VERSION
    END_VERSIONS
    """
}
