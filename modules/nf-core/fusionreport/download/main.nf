process FUSIONREPORT_DOWNLOAD {
    tag 'fusionreport'
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d9/d9d1075dc45da6b08ec99c6e8bcc83e0ab71a674e7efdc7a36e459539793fcf9/data' :
        'community.wave.seqera.io/library/fusion-report_openpyxl:6748677442b83a9a'}"

    output:
    tuple val(meta), path("fusion_report_db"), emit: fusionreport_ref
    path "versions.yml"                      , emit: versions

    script:
    meta = [id: 'fusion_report_db']
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
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
