process RRNATRANSCRIPTS {
    tag '$rrna'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(gtf)

    output:
    path("*rrna_intervals.gtf")   , emit: rrna_gtf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    template "get_rrna_transcripts.py $gtf ${prefix}_rrna_intervals.gtf"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gtf
    touch versions.yml
    """
}
