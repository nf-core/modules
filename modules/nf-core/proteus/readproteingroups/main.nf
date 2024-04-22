process PROTEUS_READPROTEINGROUPS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-503e259d7d34ce533ce66c4c8871af4ab409db6d:1e504ef71c83943061a39b6260d826b988bfa56f-0':
        'biocontainers/mulled-v2-503e259d7d34ce533ce66c4c8871af4ab409db6d:1e504ef71c83943061a39b6260d826b988bfa56f-0' }"

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
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'proteus_readproteingroups.R'
}
