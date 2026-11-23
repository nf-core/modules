// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BISMARK_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/106c5858e0e1021f5eb80619e7ec56f02e5fc857d2bdbb00570d30a5add9b5e9/data' :
        'community.wave.seqera.io/library/bismark:3.4.0--4721154cfbeba5e3' }"

    input:
    tuple val(meta), path(align_report), path(dedup_report), path(splitting_report), path(mbias)

    output:
    tuple val(meta), path("*report.{html,txt}"), emit: report
    tuple val("${task.process}"), val("bismark"), eval("bismark --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+'"), topic: versions, emit: versions_bismark

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bismark2report ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.report.txt
    touch ${prefix}.report.html
    """
}
