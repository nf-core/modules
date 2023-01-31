process AFFY_JUSTRMA {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bioconductor-affy=1.76.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-affy:1.76.0--r42hc0cfd56_2':
        'quay.io/biocontainers/bioconductor-affy:1.76.0--r42hc0cfd56_2' }"

    input:
    tuple val(meta), path(samplesheet), path(celfiles_dir)
    tuple val(meta2), path(description)

    output:
    tuple val(meta), path("*.rds"), emit: rds
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'affy_justrma.R'
}
