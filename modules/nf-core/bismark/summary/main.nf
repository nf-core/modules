process BISMARK_SUMMARY {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a65c28a667a89edfb989e4e47a01246ce75577673e02e103ab1cd30fbca84d31/data' :
        'community.wave.seqera.io/library/bismark:3.0.0--e50ebd38454e3a10' }"

    input:
    val(bam)
    path(align_report)
    path(dedup_report)
    path(splitting_report)
    path(mbias)

    output:
    path("*report.{html,txt}"), emit: summary
    tuple val("${task.process}"), val('bismark'), eval("bismark --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+'"), emit: versions_bismark, topic: versions

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
