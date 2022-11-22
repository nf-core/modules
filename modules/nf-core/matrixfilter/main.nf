process MATRIXFILTER {
    tag "$meta"
    label 'process_single'
    conda (params.enable_conda ? "r-base" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'quay.io/biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(samplesheet), path(abundance)

    output:
    tuple val(meta), path("*.filtered.tsv")             , emit: filtered
    tuple val(meta), path("R_sessionInfo.log")          , emit: session_info
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'matrixfilter.R'
}
