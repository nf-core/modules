process FUSIONREPORT_DOWNLOAD {
    tag 'fusionreport'
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a4/a41e89223340438d7cebca8af469cef9189d177bb44373b1bcb010c8e441a444/data' :
        'community.wave.seqera.io/library/fusion-report_beautifulsoup4_click_colorlog_pruned:4aeac839f8a4fc86'}"

    input:
    val(meta)

    output:
    tuple val(meta), path("${prefix}"), emit: fusionreport_ref
    tuple val("${task.process}"), val('fusion_report'), eval("fusion_report --version |& sed 's/fusion-report //'"), topic: versions, emit: versions_fusionreport


    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    fusion_report download ${args} ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/cosmic.db
    touch ${prefix}/fusiongdb2.db
    touch ${prefix}/mitelman.db
    touch ${prefix}/DB-timestamp.txt
    touch ${prefix}/fusion_report.log
    """
}
