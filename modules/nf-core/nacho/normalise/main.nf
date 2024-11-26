nextflow.enable.moduleBinaries = true

process NACHO_NORMALISE {
    tag "$sample_sheet"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/r-dplyr_r-fs_r-ggplot2_r-nacho_pruned:2d777db141950c99'

    input:
    path rcc_files, stageAs: "input/*"
    path sample_sheet

    output:
    path "*normalized_counts.tsv"          , emit: normalized_counts
    path "*normalized_counts_wo_HKnorm.tsv", emit: normalized_counts_wo_HK
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    nacho_norm.R --input_rcc_path input --input_samplesheet ${sample_sheet} $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-nacho: \$(Rscript -e "library(NACHO); cat(as.character(packageVersion('NACHO')))")
        r-dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        r-tidyr: \$(Rscript -e "library(tidyr); cat(as.character(packageVersion('tidyr')))")
        r-readr: \$(Rscript -e "library(readr); cat(as.character(packageVersion('readr')))")
        r-fs: \$(Rscript -e "library(fs); cat(as.character(packageVersion('fs')))")
        r-optparse: \$(Rscript -e "library(optparse); cat(as.character(packageVersion('optparse')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch normalized_counts.tsv
    touch normalized_counts_wo_HKnorm.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-nacho: \$(Rscript -e "library(NACHO); cat(as.character(packageVersion('NACHO')))")
        r-dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        r-tidyr: \$(Rscript -e "library(tidyr); cat(as.character(packageVersion('tidyr')))")
        r-readr: \$(Rscript -e "library(readr); cat(as.character(packageVersion('readr')))")
        r-fs: \$(Rscript -e "library(fs); cat(as.character(packageVersion('fs')))")
        r-optparse: \$(Rscript -e "library(optparse); cat(as.character(packageVersion('optparse')))")
    END_VERSIONS
    """
}
