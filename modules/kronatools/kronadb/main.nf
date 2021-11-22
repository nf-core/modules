process KRONATOOLS_KRONADB {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::krona=2.7.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/krona:2.7.1--pl526_5"
    } else {
        container "quay.io/biocontainers/krona:2.7.1--pl526_5"
    }
    input:

    output:
    path 'taxonomy/taxonomy.tab', emit: db
    path "versions.yml"         , emit: versions

    script:
    def args = task.ext.args ?: ''
    def VERSION='2.7.1'
    """
    ktUpdateTaxonomy.sh taxonomy

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: $VERSION
    END_VERSIONS
    """
}
