process GSEA_GSEA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f4fe28961396eeeaa98484cb4f2db5c79abfdf117700df132312fe5c41bff81/data':
        'community.wave.seqera.io/library/gsea:4.3.2--a7421d7504fd7c81' }"

    input:
    tuple val(meta), path(gct), path(cls), path(gene_sets)
    tuple val(reference), val(target)
    path(chip) // Optional identifier mapping file

    output:
    tuple val(meta), path("*.rpt")                             , emit: rpt
    tuple val(meta), path("*index.html")                       , emit: index_html
    tuple val(meta), path("*heat_map_corr_plot.html")          , emit: heat_map_corr_plot
    tuple val(meta), path("*gsea_report_for_${reference}.tsv") , emit: report_tsvs_ref
    tuple val(meta), path("*gsea_report_for_${reference}.html"), emit: report_htmls_ref
    tuple val(meta), path("*gsea_report_for_${target}.tsv")    , emit: report_tsvs_target
    tuple val(meta), path("*gsea_report_for_${target}.html")   , emit: report_htmls_target
    tuple val(meta), path("*ranked_gene_list*.tsv")            , emit: ranked_gene_list
    tuple val(meta), path("*gene_set_sizes.tsv")               , emit: gene_set_sizes
    tuple val(meta), path("*global_es_histogram.png")          , emit: histogram
    tuple val(meta), path("*heat_map_1.png")                   , emit: heatmap
    tuple val(meta), path("*pvalues_vs_nes_plot.png")          , emit: pvalues_vs_nes_plot
    tuple val(meta), path("*ranked_list_corr_2.png")           , emit: ranked_list_corr
    tuple val(meta), path("*butterfly_plot.png")               , emit: butterfly_plot  , optional: true
    tuple val(meta), path("gene_sets_*.tsv")                   , emit: gene_set_tsv    , optional: true
    tuple val(meta), path("gene_sets_*.html")                  , emit: gene_set_html   , optional: true
    tuple val(meta), path("gene_sets_*.png")                   , emit: gene_set_heatmap, optional: true
    tuple val(meta), path("*_snapshot*.html")                  , emit: snapshot        , optional: true
    tuple val(meta), path("*enplot*.png")                      , emit: gene_set_enplot , optional: true
    tuple val(meta), path("*gset_rnd_es_dist*.png")            , emit: gene_set_dist   , optional: true
    tuple val(meta), path("*.zip")                             , emit: archive         , optional: true
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rpt_label = prefix.replaceAll('\\.$', '') // Remove any trailing dots from prefix when passed as report label, so GSEA doesn't produce double-dotted top-level outputs
    def chip_command = chip ? "-chip $chip -collapse true" : ''
    def VERSION = '4.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    if (!(args ==~ /.*-rnd_seed.*/)) {args += " -rnd_seed 10"}

    """
    # Run GSEA

    gsea-cli GSEA \\
        -res $gct \\
        -cls ${cls}#${target}_versus_${reference} \\
        -gmx $gene_sets \\
        $chip_command \\
        -out . \\
        --rpt_label $rpt_label \\
        $args

    # Un-timestamp the outputs for path consistency
    mv ${rpt_label}.Gsea.*/* .
    timestamp=\$(cat *.rpt | grep producer_timestamp | awk '{print \$2}')

    for pattern in _\${timestamp} .\${timestamp}; do
        find . -name "*\${pattern}*" | sed "s|^\\./||" | while read -r f; do
            mv \$f \${f//\$pattern/}
        done
    done
    sed -i.bak "s/[_\\.]\$timestamp//g" *.rpt *.html && rm *.bak

    # Prefix files that currently lack it
    ls -p | grep -v / | grep -v "$prefix" | while read -r f; do
        mv \$f ${prefix}\${f}
        sed -i.bak "s/\$f/${prefix}\${f}/g" *.rpt *.html && rm *.bak
    done

    # Rename files so that they can be properly referenced by the output channels
    # Function to rename files based on the given pattern
    rename_files() {
        local pattern=\$1
        local exclude_patterns=\$2
        local extension=\$3

        # Find files matching the pattern but not matching the exclusion patterns
        find . -type f -name "\$pattern" | while read -r file; do
            # Exclude files based on the provided exclusion patterns
            if ! echo "\$file" | grep -qE "\$exclude_patterns"; then
                # Rename the file by adding the prefix "gene_sets_"
                mv "\$file" "\$(dirname "\$file")/gene_sets_\$(basename "\$file")"
            fi
        done
    }

    # Pattern and exclusion for .tsv files
    tsv_pattern="*.tsv"
    tsv_exclude="gene_set_size|gsea_report|ranked_gene_list"

    # Pattern and exclusion for .html files
    html_pattern="*.html"
    html_exclude="gsea_report|heat_map_corr_plot|index|pos_snapshot|neg_snapshot"

    # Pattern and exclusion for .png files
    png_pattern="*.png"
    png_exclude="butterfly|enplot|global_es_histogram|gset_rnd_es_dist|heat_map|pvalues_vs_nes_plot|ranked_list_corr"

    # Rename .tsv files
    rename_files "\$tsv_pattern" "\$tsv_exclude" ".tsv"

    # Rename .html files
    rename_files "\$html_pattern" "\$html_exclude" ".html"

    # Rename .png files
    rename_files "\$png_pattern" "\$png_exclude" ".png"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gsea: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '4.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.rpt
    touch ${prefix}.index.html
    touch ${prefix}.heat_map_corr_plot.html
    touch ${prefix}.gsea_report_for_${reference}.tsv
    touch ${prefix}.gsea_report_for_${reference}.html
    touch ${prefix}.gsea_report_for_${target}.tsv
    touch ${prefix}.gsea_report_for_${target}.html
    touch ${prefix}.ranked_gene_list*.tsv
    touch ${prefix}.gene_set_sizes.tsv
    touch ${prefix}.global_es_histogram.png
    touch ${prefix}.heat_map_1.png
    touch ${prefix}.pvalues_vs_nes_plot.png
    touch ${prefix}.ranked_list_corr_2.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gsea: $VERSION
    """
}
