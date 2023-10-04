process PROPR_PROPD {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::r-propr=4.2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-propr:4.2.6':
        'biocontainers/r-propr:4.2.6' }"

    input:
    tuple val(meta), path(count)
    tuple val(meta2), path(samplesheet)

    output:
    tuple val(meta), path("*.propd.tsv"), emit: results
    tuple val(meta), path("*.fdr.tsv")  , emit: fdr         , optional:true
    path "*.R_sessionInfo.log"          , emit: session_info
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'propd.R'
}
