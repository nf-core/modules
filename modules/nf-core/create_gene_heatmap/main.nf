process CREATE_GENE_HEATMAP {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bioconductor-complexheatmap_r-base_r-circlize_r-dplyr_pruned:58d1af3dbaeba617"

    input:
    tuple val(meta), path(annotated_endo_data), path(normalized_counts)
    path heatmap_genes_to_filter

    output:
    tuple val(meta), path("*gene_heatmap_mqc.png"), emit: gene_heatmap
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    compute_gene_heatmap.R $annotated_endo_data $normalized_counts $heatmap_genes_to_filter $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        r-rlang: \$(Rscript -e "library(rlang); cat(as.character(packageVersion('rlang')))")
        bioconductor-ComplexHeatmap: \$(Rscript -e "library(ComplexHeatmap); cat(as.character(packageVersion('ComplexHeatmap')))")
        r-circlize: \$(Rscript -e "library(circlize); cat(as.character(packageVersion('circlize')))")
        r-yaml: \$(Rscript -e "library(yaml); cat(as.character(packageVersion('yaml')))")
        r-fs: \$(Rscript -e "library(fs); cat(as.character(packageVersion('fs')))")
    END_VERSIONS
    """

    stub:
    """
    touch gene_heatmap_mqc.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        r-rlang: \$(Rscript -e "library(rlang); cat(as.character(packageVersion('rlang')))")
        bioconductor-ComplexHeatmap: \$(Rscript -e "library(ComplexHeatmap); cat(as.character(packageVersion('ComplexHeatmap')))")
        r-circlize: \$(Rscript -e "library(circlize); cat(as.character(packageVersion('circlize')))")
        r-yaml: \$(Rscript -e "library(yaml); cat(as.character(packageVersion('yaml')))")
        r-fs: \$(Rscript -e "library(fs); cat(as.character(packageVersion('fs')))")
    END_VERSIONS
    """

}

