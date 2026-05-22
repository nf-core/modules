process DOTSEQ_DOTSEQ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/12/12667d472e9ae0f1602041dc018ba6bde294e6190e67999d71b65e7a2df7ea1f/data' :
        'community.wave.seqera.io/library/bioconductor-dotseq_r-dplyr_r-eulerr_r-ggplot2_pruned:6c8a9ebdec36c958' }"

    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target)
    tuple val(meta2), path(samplesheet), path(counts), path(annotation)

    output:
    tuple val(meta), path("*.translation.dotseq.results.tsv")    , emit: translation
    tuple val(meta), path("*.dou.dotseq.results.tsv")            , emit: dou
    tuple val(meta), path("*.dou_strategy.dotseq.results.tsv")   , emit: dou_strategy , optional: true
    tuple val(meta), path("*.dte_strategy.dotseq.results.tsv")   , emit: dte_strategy , optional: true
    tuple val(meta), path("*.volcano.png")                       , emit: volcano_plot  , optional: true
    tuple val(meta), path("*.composite.png")                     , emit: composite_plot, optional: true
    tuple val(meta), path("*.venn.png")                          , emit: venn_plot     , optional: true
    tuple val(meta), path("*.heatmap.png")                       , emit: heatmap_plot  , optional: true
    tuple val(meta), path("*.interaction_p_distribution.png")    , emit: interaction_p_distribution_plot, optional: true
    tuple val(meta), path("*.DOTSeqDataSets.rds")                , emit: rdata
    tuple val(meta), path("*.R_sessionInfo.log")                 , emit: session_info
    path "versions.yml"                                          , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'dotseq.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.translation.dotseq.results.tsv
    touch ${prefix}.dou.dotseq.results.tsv
    touch ${prefix}.dou_strategy.dotseq.results.tsv
    touch ${prefix}.dte_strategy.dotseq.results.tsv
    touch ${prefix}.volcano.png
    touch ${prefix}.composite.png
    touch ${prefix}.venn.png
    touch ${prefix}.heatmap.png
    touch ${prefix}.interaction_p_distribution.png
    touch ${prefix}.DOTSeqDataSets.rds
    touch ${prefix}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-dotseq: \$(Rscript -e "cat(as.character(packageVersion('DOTSeq')))")
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        r-readr: \$(Rscript -e "cat(as.character(packageVersion('readr')))")
        r-dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
