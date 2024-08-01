process SEURAT_MULTISEQ {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-seurat:3.0.2--r36h0357c0b_1':
        'quay.io/biocontainers/r-seurat:3.0.2--r36h0357c0b_0' }"

    input:
    tuple val(meta), path(hto_matrix),path(rna_matrix)
    output:
    tuple val(meta), path("*.csv")                                  ,emit: assignment
    tuple val(meta), path("RidgePlot.jpeg")         ,optional:true  ,emit: ridgePlot
    tuple val(meta), path("FeatureScatter.jpeg")    ,optional:true  ,emit: featureScatter
    tuple val(meta), path("ViolinPlot.jpeg")        ,optional:true  ,emit: violingPlot
    tuple val(meta), path("tSNE.jpeg")              ,optional:true  ,emit: tSNEPlot
    path "versions.yml"                             .optional:true  ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    template 'MULTIseq.R'
}
