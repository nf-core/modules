process BISMARK_SUMMARY {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/38e61d14ccaed82f60c967132963eb467d0fa4bccb7a21404c49b4f377735f03/data' :
        'community.wave.seqera.io/library/bismark:0.25.1--1f50935de5d79c47' }"

    input:
    val(bam)
    path(align_report)
    path(dedup_report)
    path(splitting_report)
    path(mbias)

    output:
    path("*report.{html,txt}"), emit: summary
    tuple val("${task.process}"), val('bismark'), eval('bismark --version 2>&1 | sed -n "s/^.*Bismark Version: v//p" | tr -d "[:space:]"'), emit: versions_bismark, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bismark2summary ${bam.join(' ')}
    """

    stub:
    """
    touch bismark_summary_report.txt
    touch bismark_summary_report.html
    """
}
