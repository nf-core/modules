process NACHO_QC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/r-dplyr_r-fs_r-ggplot2_r-nacho_pruned:033bc017f5f36b6d'

    input:
    tuple val(meta) , path(rcc_files, stageAs: "input/*")
    tuple val(meta2), path(sample_sheet)

    output:
    tuple val(meta), path("*.html")   , emit: nacho_qc_reports
    tuple val(meta), path("*_mqc.png"), emit: nacho_qc_png
    tuple val(meta), path("*_mqc.txt"), emit: nacho_qc_txt
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    nacho_qc.R \\
        --input_rcc_path input \\
        --input_samplesheet ${sample_sheet}

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
    """
    touch qc.html
    touch qc_with_outliers.html
    touch AVG_vs_BD_mqc.png
    touch AVG_vs_MED_mqc.png
    touch BD_mqc.png
    touch FOV_mqc.png
    touch HKF_mqc.png
    touch HK_mqc.png
    touch LOD_mqc.png
    touch Neg_mqc.png
    touch PCA1_vs_PCA2_mqc.png
    touch PCAi_mqc.png
    touch PCA_mqc.png
    touch plot_normf_mqc.png
    touch Posctrl_linearity_mqc.png
    touch POSF_vs_NEGF_mqc.png
    touch Pos_mqc.png
    touch Pos_vs_neg_mqc.png
    touch normalized_qc_mqc.txt
    touch hk_detected_mqc.txt

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
