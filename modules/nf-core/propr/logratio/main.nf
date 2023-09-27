process PROPR_LOGRATIO {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::r-propr=4.2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-propr:4.2.6':
        'biocontainers/r-propr:4.2.6' }"

    input:
    tuple val(meta), path(count), val(transformation), val(reference), val(alpha)

    output:
    tuple val(meta), path("*.logratio.tsv"), emit: logratio
    path "*.R_sessionInfo.log", emit: session_info
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'logratio.R'
}