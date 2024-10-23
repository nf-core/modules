process PROPR_PROPR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-propr:5.0.3':
        'biocontainers/r-propr:5.0.3' }"

    input:
    tuple val(meta), path(count)

    output:
    tuple val(meta), path("*.propr.rds"), emit: propr
    tuple val(meta), path("*.propr.tsv"), emit: matrix
    tuple val(meta), path("*.fdr.tsv"),   emit: fdr         , optional:true
    tuple val(meta), path("*.adj.csv"),   emit: adj         , optional:true
    path "*.warnings.log",                emit: warnings
    path "*.R_sessionInfo.log",           emit: session_info
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'propr.R'
}
