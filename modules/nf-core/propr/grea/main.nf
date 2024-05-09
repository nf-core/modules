process PROPR_GREA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-propr:5.0.4':
        'biocontainers/r-propr:5.0.4' }"

    input:
    tuple val(meta), path(adj)
    tuple val(meta2), path(gmt)

    output:
    tuple val(meta), path("*.go.tsv"),  emit: enrichedGO
    path "versions.yml",                emit: versions
    path "*.R_sessionInfo.log",         emit: session_info

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'grea.R'
}
