def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/fcsgx/rungx

Reason:
This module is now renamed as FCSGX_RUNGX and as been updated to the latest version
"""

process FCS_FCSGX {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(assembly)
    path gxdb

    output:
    tuple val(meta), path("out/*.fcs_gx_report.txt"), emit: fcs_gx_report
    tuple val(meta), path("out/*.taxonomy.rpt")     , emit: taxonomy_report
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    assert false: deprecation_message

    stub:
    assert false: deprecation_message
}
