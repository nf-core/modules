process NACHO_QC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ad/ad421367a5d71eb73738675a68b5677e283686a8b0a6d5e5530f9ec203aadb30/data' :
        'community.wave.seqera.io/library/r-base_r-dplyr_r-fs_r-ggplot2_pruned:bcd6d91c836e9200' }"

    input:
    tuple val(meta) , path(rcc_files, stageAs: "input/*")
    tuple val(meta2), path(sample_sheet)

    output:
    tuple val(meta), path("*.html")   , emit: nacho_qc_reports
    tuple val(meta), path("*_mqc.png"), emit: nacho_qc_png
    tuple val(meta), path("*_mqc.txt"), emit: nacho_qc_txt
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    template 'nacho_qc.R'

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

    Rscript -e "nfcore.utils::process_end(
        packages = list(
            'r-nacho' = 'NACHO',
            'r-dplyr' = 'dplyr',
            'r-ggplot2' = 'ggplot2',
            'r-tidyr' = 'tidyr',
            'r-readr' = 'readr',
            'r-fs' = 'fs'
        ),
        task_name = '${task.process}',
        versions_path = 'versions.yml',
        log_path = 'R_sessionInfo.log'
    )"
    """
}
