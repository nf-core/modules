process KRONA_KRONADB {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krona:2.7.1--pl526_5' :
        'quay.io/biocontainers/krona:2.7.1--pl526_5' }"

    output:
    path 'taxonomy/taxonomy.tab', emit: db
    tuple val("${task.process}"), val('krona'), eval("ktImportTaxonomy | grep -Po \"(?<=KronaTools )[0-9.]+\""), emit: versions_krona, topic: versions

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
    """
    ktUpdateTaxonomy.sh \\
        $args \\
        taxonomy/

    """

    stub:
    """
    mkdir taxonomy

    touch taxonomy/taxonomy.tab

    """
}
