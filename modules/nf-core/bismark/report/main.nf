process BISMARK_REPORT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bddea334e6ccbce005ce540214747acf822b040185d2198220dcfbb4b258c331/data'
        : 'community.wave.seqera.io/library/bismark:3.1.0--9557d6ab108a83e4'}"

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
