process CUSTOM_MATRIXFILTER {
    tag "$meta"
    label 'process_single'
    conda (params.enable_conda ? "r-base" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'quay.io/biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(abundance)
    tuple val(samplesheet_meta), path(samplesheet)

    output:
    tuple val(meta), path("*.filtered.tsv")             , emit: filtered
    tuple val(meta), path("R_sessionInfo.log")          , emit: session_info
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Note: params are passed at line 100 of the template like:
    //
    // opt <- parse_args('$task.ext.args', opt)
    //
    // (new variables defined here don't seem to be available in templates, so
    // we have to access $task directly)
    template 'matrixfilter.R'
}
