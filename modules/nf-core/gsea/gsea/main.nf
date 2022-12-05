process GSEA_GSEA {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::gsea=4.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gsea:4.3.2--hdfd78af_0':
        'quay.io/biocontainers/gsea:4.3.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(gct), path(cls), path(gene_sets)

    output:
    tuple val(meta), path("*/*.rpt")                                                                                                        , emit: param_log
    tuple val(meta), path("*/index.html")                                                                                                   , emit: index
    tuple val(meta), path("*/heat_map_corr_plot.html")                                                                                      , emit: heat_map_corr_plot
    tuple val(meta), path("*/gsea_report_for_${meta.reference}*.tsv"),  path("*/gsea_report_for_${meta.target}*.tsv")                       , emit: report_tsvs
    tuple val(meta), path("*/gsea_report_for_${meta.reference}*.html"),  path("*/gsea_report_for_${meta.target}*.html")                     , emit: report_htmls
    tuple val(meta), path("*/ranked_gene_list*.tsv")                                                                                        , emit: ranked_gene_list
    tuple val(meta), path("*/gene_set_sizes.tsv")                                                                                           , emit: gene_set_sizes
    tuple val(meta), path("*/butterfly_plot.png")                                                                                           , emit: butterfly_plot
    tuple val(meta), path("*/global_es_histogram.png")                                                                                      , emit: histogram
    tuple val(meta), path("*/heat_map_1.png")                                                                                               , emit: heatmap
    tuple val(meta), path("*/pvalues_vs_nes_plot.png")                                                                                      , emit: pvalues_vs_nes_plot
    tuple val(meta), path("*/ranked_list_corr_2.png")                                                                                       , emit: ranked_list_corr
    tuple val(meta), path("*/[!gene_set_size|gsea_report|ranked_gene_list]*.tsv")                                                           , emit: gene_set_tsv, optional: true
    tuple val(meta), path("*/[!gsea_report|heat_map_corr_plot|index|pos_snapshot|neg_snapshot]*.html")                                      , emit: gene_set_html, optional: true
    tuple val(meta), path("*/[!butterfly|enplot|global_es_histogram|gset_rnd_es_dist|heat_map|pvalues_vs_nes_plot|ranked_list_corr]*.png")  , emit: gene_set_heatmap, optional: true
    tuple val(meta), path("*/*_snapshot*.html")                                                                                             , emit: snapshot, optional: true
    tuple val(meta), path("*/enplot*.png")                                                                                                  , emit: gene_set_enplot, optional: true
    tuple val(meta), path("*/gset_rnd_es_dist*.png")                                                                                        , emit: gene_set_dist, optional: true
    tuple val(meta), path("*/*.zip")                                                                                                        , emit: archive, optional: true
    path "versions.yml"                                                                                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Run GSEA
    gsea-cli GSEA \
        -res $gct \
        -cls ${cls}#${meta.target}_versus_${meta.reference} \
        -gmx $gene_sets \
        -out \$(pwd) \
        --rpt_label $prefix \
        $args

    # Move things out of the timestamped folder
    mv ${prefix}.Gsea.* ${prefix}

    # Other files have the annoying timestamp prefix too, but may be referred
    # to from other files (e.g. the index.html). Symlinked here to provide
    # consistent paths for the nf-core CI

    pushd $prefix
    ln -s ${prefix}.Gsea.*.rpt ${prefix}.Gsea.rpt
    (ls *.zip 2>/dev/null || true) | while read -r l; do
        ln -s \$l ${prefix}.Gsea.rpt.zip
    done

    ln -s \$(ls ranked_gene_list_${meta.target}_versus_${meta.reference}_*.tsv) ranked_gene_list_${meta.target}_versus_${meta.reference}.tsv
    ln -s \$(ls gsea_report_for_${meta.reference}_*.html) gsea_report_for_${meta.reference}.html
    ln -s \$(ls gsea_report_for_${meta.reference}_*.tsv) gsea_report_for_${meta.reference}.tsv
    ln -s \$(ls gsea_report_for_${meta.target}_*.html) gsea_report_for_${meta.target}.html
    ln -s \$(ls gsea_report_for_${meta.target}_*.tsv) gsea_report_for_${meta.target}.tsv
    popd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gsea: 4.3.2
    END_VERSIONS
    """
}
