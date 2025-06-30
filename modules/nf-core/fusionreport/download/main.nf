process FUSIONREPORT_DOWNLOAD {
    tag 'fusionreport'
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/fusion-report_beautifulsoup4_click_colorlog_pruned:353260d00d0da0ec' :
        'community.wave.seqera.io/library/fusion-report_beautifulsoup4_click_colorlog_pruned:78488bd99166aa9a'}"

    output:
    tuple val(meta), path("fusion_report_db"), emit: fusionreport_ref
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    meta = [id: 'fusion_report_db']
    def args = task.ext.args ?: ''
    """
    fusion_report download $args ./
    mkdir fusion_report_db
    mv *.txt *.log *.db fusion_report_db/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusion_report: \$(fusion_report --version | sed 's/fusion-report //')
    END_VERSIONS
    """

    stub:
    meta = [id: 'fusion_report_db']
    """
    mkdir fusion_report_db
    touch fusion_report_db/cosmic.db
    touch fusion_report_db/fusiongdb2.db
    touch fusion_report_db/mitelman.db
    touch fusion_report_db/DB-timestamp.txt
    touch fusion_report_db/fusion_report.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusion_report: \$(fusion_report --version | sed 's/fusion-report //')
    END_VERSIONS
    """
}
