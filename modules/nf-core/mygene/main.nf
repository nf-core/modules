process MYGENE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mygene:3.2.2--pyh5e36f6f_0':
        'biocontainers/mygene:3.2.2--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(gene_list)

    output:
    tuple val(meta), path("*.gmt"), emit: gmt
    tuple val(meta), path("*.tsv"), emit: tsv     , optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "mygene.py"
}
