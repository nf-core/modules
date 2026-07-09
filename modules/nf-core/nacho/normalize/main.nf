process NACHO_NORMALIZE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9e04cca7290d2e7b649e4523fa7aa72660d201908c8003c5690546fe8dac243d/data' :
        'community.wave.seqera.io/library/r-dplyr_r-fs_r-ggplot2_r-nacho_pruned:48ffc9bbb33907fc' }"

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
