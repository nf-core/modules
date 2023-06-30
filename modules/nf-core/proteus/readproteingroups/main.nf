process PROTEUS_READPROTEINGROUPS {
    tag "$meta"
    label 'process_single'

    conda "r-base=4.2.1 r-proteus-bartongroup=0.2.16 bioconductor-limma=3.54.0 r-plotly=4.10.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4e01206f2c47f56077f04e5d2d7b312f50513a1e:92abccefbeb09795ad6a93553b62a6ad3daaea48-0':
        'quay.io/biocontainers/mulled-v2-4e01206f2c47f56077f04e5d2d7b312f50513a1e:92abccefbeb09795ad6a93553b62a6ad3daaea48-0' }"

    input:
    tuple val(meta), path(samplesheet), path(intensities)
    tuple val(meta2), val(contrast_variable)


    output:
    tuple val(meta), path("*dendrogram.png")                    , emit: dendro_plot
    tuple val(meta), path("*mean_variance_relationship.png")    , emit: mean_var_plot
    tuple val(meta), path("*normalised_distributions.png")      , emit: raw_dist_plot
    tuple val(meta), path("*normalised_distributions.png")      , emit: norm_dist_plot
    tuple val(meta), path("*raw_proteingroups.rds")             , emit: raw_rdata
    tuple val(meta), path("*raw_proteingroups_tab.tsv")         , emit: raw_tab
    tuple val(meta), path("*normalised_proteingroups.rds")      , emit: norm_rdata
    tuple val(meta), path("*normalised_proteingroups_tab.tsv")  , emit: norm_tab
    tuple val(meta), path("*R_sessionInfo.log")                 , emit: session_info
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'proteus_readproteingroups.R'
}
