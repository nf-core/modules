process ISOSEQ3_TAG {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoseq3:4.0.0--h9ee0642_0':
        'quay.io/biocontainers/isoseq3:4.0.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam)
    val design

    output:
    tuple val(meta), path("*.flt.bam")                  , emit: bam
    tuple val(meta), path("*.flt.bam.pbi")              , emit: pbi
    tuple val("${task.process}"), val('isoseq3'), eval("isoseq tag --version | sed -n '1s/isoseq tag \\([0-9.]*\\).*/\\1/p'"), topic: versions, emit: versions_isoseq3

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
    WARNING: This module has been deprecated. Please use nf-core/modules/isoseq/tag

    Reason:
    The isoseq3 conda package is frozen at 4.0.0 and now only wraps the actively
    maintained isoseq package. isoseq/tag provides the same functionality under
    the current isoseq release.
    """
    assert false: deprecation_message

    stub:
    def deprecation_message = """
    WARNING: This module has been deprecated. Please use nf-core/modules/isoseq/tag

    Reason:
    The isoseq3 conda package is frozen at 4.0.0 and now only wraps the actively
    maintained isoseq package. isoseq/tag provides the same functionality under
    the current isoseq release.
    """
    assert false: deprecation_message
}
