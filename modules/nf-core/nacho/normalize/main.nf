process NACHO_NORMALIZE {
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
    tuple val(meta), path("${prefix}.tsv")          , emit: normalized_counts
    tuple val(meta), path("${prefix}_wo_HKnorm.tsv"), emit: normalized_counts_wo_HK
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo ${args}
    """

    template 'nacho_norm.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    args = task.ext.args ?: ''

    """
    touch ${prefix}.tsv
    touch ${prefix}_wo_HKnorm.tsv

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
        log_path = '${prefix}.R_sessionInfo.log'
    )"
    """
}
