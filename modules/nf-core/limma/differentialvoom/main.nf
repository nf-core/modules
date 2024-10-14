process LIMMA_DIFFERENTIALVOOM {
    tag "$meta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-edger_bioconductor-ihw_bioconductor-limma_r-dplyr_r-readr:7fc48564d286c1c6' :
        'community.wave.seqera.io/library/bioconductor-edger_bioconductor-ihw_bioconductor-limma_r-dplyr_r-readr:edea0f9fbaeba3a0' }"

    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target)
    tuple val(meta2), path(samplesheet), path(intensities)

    output:
    tuple val(meta), path("*.limma.results.tsv")          , emit: results
    tuple val(meta), path("*.limma.mean_difference.png")  , emit: md_plot
    tuple val(meta), path("*.MArrayLM.limma.rds")         , emit: rdata
    tuple val(meta), path("*.normalised_counts.tsv")      , emit: normalised_counts, optional: true
    tuple val(meta), path("*.limma.model.txt")            , emit: model
    tuple val(meta), path("*.R_sessionInfo.log")          , emit: session_info
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'limma_de.R'

    stub:
    prefix              = task.ext.prefix   ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    library(limma)

    a <- file("${prefix}.limma.results.tsv", "w")
    close(a)
    a <- file("${prefix}.limma.mean_difference.png", "w")
    close(a)
    a <- file("${prefix}.MArrayLM.limma.rds", "w")
    close(a)
    a <- file("${prefix}.normalised_counts.tsv", "w")
    close(a)
    a <- file("${prefix}.limma.model.txt", "w")
    close(a)
    a <- file("${prefix}.R_sessionInfo.log", "w")
    close(a)

    ## VERSIONS FILE
    r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
    limma.version <- as.character(packageVersion('limma'))

    writeLines(
        c(
            '"task.process":',
            paste('    r-base:', r.version),
            paste('    bioconductor-limma:', limma.version)
        ),
        'versions.yml'
    )
    """
}
