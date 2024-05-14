process LIMMA_DIFFERENTIAL {
    tag "$meta"
    label 'process_medium'

    conda "bioconda::bioconductor-limma=3.54.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-limma:3.54.0--r42hc0cfd56_0' :
        'biocontainers/bioconductor-limma:3.54.0--r42hc0cfd56_0' }"

    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target)
    tuple val(meta2), path(samplesheet), path(intensities)

    output:
    tuple val(meta), path("*.limma.results.tsv")          , emit: results
    tuple val(meta), path("*.limma.mean_difference.png")  , emit: md_plot
    tuple val(meta), path("*.MArrayLM.limma.rds")         , emit: rdata
    tuple val(meta), path("*.R_sessionInfo.log")          , emit: session_info
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'limma_de.R'
}
