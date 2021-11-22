process KRONATOOLS_KTIMPORTTAXONOMY {
    tag "${meta.id}"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::krona=2.8" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/krona:2.8--pl5262hdfd78af_2"
    } else {
        container "quay.io/biocontainers/krona:2.8--pl5262hdfd78af_2"
    }

    input:
    tuple val(meta), path(report)
    path  "taxonomy/taxonomy.tab"

    output:
    tuple val(meta), path ('*.html'), emit: html
    path "versions.yml"             , emit: versions

    script:
    def VERSION='2.8'
    """
    ktImportTaxonomy "$report" -tax taxonomy

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: $VERSION
    END_VERSIONS
    """
}
