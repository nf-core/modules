// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BISMARK_SUMMARY {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/65/65d4877a6dc7b42a4fb5086d4620cdfdba48a5e2407e81f1de27357743b4aaf1/data' :
        'community.wave.seqera.io/library/bismark:3.2.0--3978d640d774e7f9' }"

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
