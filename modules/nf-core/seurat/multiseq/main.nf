process SEURAT_MULTISEQ {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-seurat:3.0.2--r36h0357c0b_1':
        'quay.io/biocontainers/r-seurat:3.0.2--r36h0357c0b_0' }"

    input:
    tuple val(meta), path(hto_matrix)
    tuple val(meta), path(rna_matrix)
    val produce_plots
    output:
    tuple val(meta), path("*.csv")                  ,emit: assignment
    path "versions.yml"                             ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    template 'MULTIseq.R'
}
