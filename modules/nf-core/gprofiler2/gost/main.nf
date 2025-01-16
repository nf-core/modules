process GPROFILER2_GOST {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e4/e4b0e10a72db4ad519c128c4e3cef6e10bc1a83440af31f105ab389a5532589a/data':
        'community.wave.seqera.io/library/r-ggplot2_r-gprofiler2:fab855ea9f680400' }"

    input:
    tuple val(meta), path(de_file)
    path(gmt_file)
    path(background_file)

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

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.gprofiler2.all_enriched_pathways.tsv
    touch ${prefix}.gprofiler2.gost_results.rds
    touch ${prefix}.gprofiler2.gostplot.png
    touch ${prefix}.gprofiler2.gostplot.html
    touch ${prefix}.gprofiler2.*.sub_enriched_pathways.tsv
    touch ${prefix}.gprofiler2.*.sub_enriched_pathways.png
    touch ${prefix}.ENSG_filtered.gmt
    touch ${prefix}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        r-gprofiler2: \$(Rscript -e "library(gprofiler2); cat(as.character(packageVersion('gprofiler2')))")
    END_VERSIONS
    """
}
