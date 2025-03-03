process KRONA_KRONADB {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krona:2.7.1--pl526_5' :
        'biocontainers/krona:2.7.1--pl526_5' }"

    output:
    path 'taxonomy/taxonomy.tab', emit: db
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
    WARNING: This module has been deprecated. Please use nf-core/modules/krona/ktupdatetaxonomy instead.

    Reason:
    This module has been superseded by the krona/ktupdatetaxonomy module.
    """
    assert false: deprecation_message

    def args = task.ext.args ?: ''
    def VERSION = '2.7.1' // Version information not provided by tool on CLI
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
