def VERSION='2.7.1' // Version information not provided by tool on CLI

process KRONA_KRONADB {
    label 'process_single'

    conda "bioconda::krona=2.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krona:2.7.1--pl526_5' :
        'quay.io/biocontainers/krona:2.7.1--pl526_5' }"

    output:
    path 'taxonomy/taxonomy.tab', emit: db
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    args += task.ext.custom_args ? ' ' + task.ext.custom_args : ''
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
