process AFFY_JUSTRMA {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bioconductor-affy=1.78.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-affy:1.78.0--r43ha9d7317_1':
        'biocontainers/bioconductor-affy:1.78.0--r43ha9d7317_1' }"

    input:
    tuple val(meta), path(samplesheet), path(celfiles_dir)
    tuple val(meta2), path(description)

    output:
    tuple val(meta), path("*.rds")             , emit: rds
    tuple val(meta), path("*matrix.tsv")       , emit: expression
    tuple val(meta), path("*.annotation.tsv")  , emit: annotation, optional: true
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'affy_justrma.R'
}
