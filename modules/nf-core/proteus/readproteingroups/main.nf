process PROTEUS_READPROTEINGROUPS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d0/d0008e2c0054e7cb33f3acef6ea417d157945a6f5774f6b7bdc1ef19cf95b249/data':
        'community.wave.seqera.io/library/bioconductor-limma_r-proteus-bartongroup_r-base_r-ggplot2_r-plotly:855c86e2407bbedb' }"

    input:
    tuple val(meta), path(samplesheet), path(intensities)

    output:
    tuple val(meta), path("*dendrogram.png")                    , emit: dendro_plot
    tuple val(meta), path("*mean_variance_relationship.png")    , emit: mean_var_plot
    tuple val(meta), path("*raw_distributions.png")             , emit: raw_dist_plot
    tuple val(meta), path("*normalized_distributions.png")      , emit: norm_dist_plot
    tuple val(meta), path("*raw_proteingroups.rds")             , emit: raw_rdata
    tuple val(meta), path("*normalized_proteingroups.rds")      , emit: norm_rdata
    tuple val(meta), path("*raw_proteingroups_tab.tsv")         , emit: raw_tab
    tuple val(meta), path("*normalized_proteingroups_tab.tsv")  , emit: norm_tab
    tuple val(meta), path("*R_sessionInfo.log")                 , emit: session_info
    path "versions.yml"                                         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'proteus_readproteingroups.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_dendrogram.png
    touch ${prefix}_mean_variance_relationship.png
    touch ${prefix}_raw_distributions.png
    touch ${prefix}_normalized_distributions.png
    touch ${prefix}_raw_proteingroups.rds
    touch ${prefix}_normalized_proteingroups.rds
    touch ${prefix}_raw_proteingroups_tab.tsv
    touch ${prefix}_normalized_proteingroups_tab.tsv
    touch ${prefix}_R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-proteus-bartongroup: \$(Rscript -e "library(r-proteus-bartongroup); cat(as.character(packageVersion('r-proteus-bartongroup')))")
        r-plotly: \$(Rscript -e "library(r-plotly); cat(as.character(packageVersion('r-plotly')))")
        bioconductor-limma: \$(Rscript -e "library(bioconductor-limma); cat(as.character(packageVersion('bioconductor-limma')))")
    END_VERSIONS
    """
}
