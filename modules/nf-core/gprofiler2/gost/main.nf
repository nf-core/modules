process GPROFILER2_GOST {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-3712554873398d849d0d11b22440f41febbc4ede:aa19bb8afc0ec6456a4f3cd650f7577c3bbdd4f3-0':
        'biocontainers/mulled-v2-3712554873398d849d0d11b22440f41febbc4ede:aa19bb8afc0ec6456a4f3cd650f7577c3bbdd4f3-0' }"

    input:
    tuple val(meta), path(de_file)
    tuple val(meta2), path(gmt_file)
    tuple val(meta3), path(background_file)

    output:
    tuple val(meta), path("*.gprofiler2.all_enriched_pathways.tsv")     , emit: all_enrich
    tuple val(meta), path("*.gprofiler2.gost_results.rds")              , emit: rds         , optional: true
    tuple val(meta), path("*.gprofiler2.gostplot.png")                  , emit: plot_png    , optional: true
    tuple val(meta), path("*.gprofiler2.gostplot.html")                 , emit: plot_html   , optional: true
    tuple val(meta), path("*.gprofiler2.*.sub_enriched_pathways.tsv")   , emit: sub_enrich  , optional: true
    tuple val(meta), path("*.gprofiler2.*.sub_enriched_pathways.png")   , emit: sub_plot    , optional: true
    tuple val(meta), path("*ENSG_filtered.gmt")                         , emit: filtered_gmt, optional: true
    tuple val(meta), path("*R_sessionInfo.log")                         , emit: session_info
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'gprofiler2_gost.R'
}
