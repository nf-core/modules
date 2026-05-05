process FUSIONREPORT_DOWNLOAD {
    tag 'fusionreport'
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed32f46746a5b33a1b5d597cfe2f62d3b1cfa54638f57cadc5d5158d6a4959d2/data' :
        'community.wave.seqera.io/library/fusion-report_beautifulsoup4_click_colorlog_pruned:78488bd99166aa9a'}"

    input:
    val(meta)

    output:
    tuple val(meta), path("fusion_report_db"), emit: fusionreport_ref
    tuple val("${task.process}"), val('fusion_report'), eval("fusion_report --version |& sed 's/fusion-report //'"), topic: versions, emit: versions_fusionreport

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    fusion_report download ${args} ./
    mkdir fusion_report_db
    mv *.txt *.log *.db fusion_report_db/
    """

    stub:
    """
    mkdir fusion_report_db
    touch fusion_report_db/cosmic.db
    touch fusion_report_db/fusiongdb2.db
    touch fusion_report_db/mitelman.db
    touch fusion_report_db/DB-timestamp.txt
    touch fusion_report_db/fusion_report.log
    """
}
