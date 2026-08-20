process MODKIT_LOCALIZE_PLOT {
    tag "$meta.id"
    label 'process_single'

    // R packages (r-base, r-data.table, r-ggplot2, r-zoo) via Seqera Wave.
    // URI built from environment.yml — freeze pending Docker Hub credentials setup.
    // Track 2: replace with stable community.wave.seqera.io URI once frozen.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://wave.seqera.io/wt/ea44ebb9e07a/wave/build:r-base-4.4.3_r-data.table-1.16.4_r-ggplot2-3.5.1_r-zoo-1.8_12--1b699ccd4840da3d' :
        'wave.seqera.io/wt/ea44ebb9e07a/wave/build:r-base-4.4.3_r-data.table-1.16.4_r-ggplot2-3.5.1_r-zoo-1.8_12--1b699ccd4840da3d' }"

    input:
    tuple val(meta), path(tsvs)
    tuple val(meta2), path(samplesheet)

    output:
    tuple val(meta), path("figures/png/*.png") , emit: png  , optional: true
    tuple val(meta), path("figures/pdf/*.pdf") , emit: pdf  , optional: true
    tuple val(meta), path("figures/svg/*.svg") , emit: svg  , optional: true
    tuple val(meta), path("*_combined.tsv")    , emit: combined_tsv , optional: true
    tuple val(meta), path("*_summary.tsv")     , emit: summary_tsv , optional: true
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'plot_localize_composite.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    dir.create('figures/png', recursive = TRUE)
    dir.create('figures/pdf', recursive = TRUE)
    dir.create('figures/svg', recursive = TRUE)
    file.create('figures/png/${prefix}_composite_lines_loess.png')
    file.create('figures/pdf/${prefix}_composite_lines_loess.pdf')
    file.create('figures/svg/${prefix}_composite_lines_loess.svg')
    file.create('${prefix}_combined.tsv')
    file.create('${prefix}_summary.tsv')

    r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', r.version)
        ),
        'versions.yml'
    )
    """
}
