process AFFY_JUSTRMA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
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
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'affy_justrma.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_eset.rds
    touch ${prefix}_matrix.tsv
    touch R_sessionInfo.log
    touch versions.yml
    """
}
