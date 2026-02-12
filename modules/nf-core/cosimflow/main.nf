process COSIMFLOW {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ad9dd5f398966bf899ae05f8e7c54d0d94018089:00bd7e543b7135f299839496732948679f2e3073-0' :
        'biocontainers/mulled-v2-ad9dd5f398966bf899ae05f8e7c54d0d94018089:00bd7e543b7135f299839496732948679f2e3073-0' }"

    input:
    tuple val(meta), path(expression_matrix)

    output:
    tuple val(meta), path("*_matrix.csv") , emit: matrix
    tuple val(meta), path("*_heatmap.png"), emit: heatmap
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages({
        library(optparse)
        library(readr)
        library(readxl)
        library(dplyr)
        library(lsa)
        library(pheatmap)
    })

    option_list = list(
        make_option(c("-i","--input"), type="character", help="Input file"),
        make_option(c("-o","--out_prefix"), type="character", default="cosine"),
        make_option(c("-m","--method"), type="character", default="cosine"),
        make_option(c("-g","--min_gene_mean"), type="double", default=0.0),
        make_option(c("-s","--sample_cols"), type="character", default="auto")
    )

    # AQUÍ PASAMOS LOS ARGUMENTOS DE NEXTFLOW DIRECTAMENTE
    args_list <- c("--input", "$expression_matrix", "--out_prefix", "$prefix")
    
    # Añadimos los argumentos extra si existen
    extra_args <- "$args"
    if (nchar(extra_args) > 0) {
        args_list <- c(args_list, unlist(strsplit(extra_args, "\\\\s+")))
    }

    opt <- parse_args(OptionParser(option_list=option_list), args=args_list)

    # ---- Load table ----
    if (grepl("\\\\.xlsx?\$", opt\$input, ignore.case=TRUE)) {
        df <- as.data.frame(read_excel(opt\$input))
    } else {
        df <- as.data.frame(read_csv(opt\$input, col_types = cols(.default = col_guess())))
    }

    # ---- Set rownames ----
    first_col <- names(df)[1]
    if (tolower(first_col) %in% c("gene","id", "gene_id")) {
        if (is.character(df[[first_col]]) || is.factor(df[[first_col]])) {
            rownames(df) <- df[[first_col]]
            df[[first_col]] <- NULL
        }
    }

    # ---- Select sample columns ----
    if (opt\$sample_cols != "auto") {
        cols <- trimws(unlist(strsplit(opt\$sample_cols,",")))
        mat <- as.matrix(df[, cols, drop=FALSE])
    } else {
        num_cols <- sapply(df, is.numeric)
        mat <- as.matrix(df[, num_cols, drop=FALSE])
    }

    # ---- Filter ----
    if (opt\$min_gene_mean > 0) {
        keep <- rowMeans(mat, na.rm=TRUE) >= opt\$min_gene_mean
        mat <- mat[keep, , drop=FALSE]
    }

    # ---- Compute ----
    method <- tolower(opt\$method)
    if (method == "cosine") {
        sim <- lsa::cosine(mat)
    } else {
        sim <- cor(mat, method = method, use = "pairwise.complete.obs")
    }

    colnames(sim) <- colnames(mat)
    rownames(sim) <- colnames(mat)

    # ---- Save ----
    write.csv(sim, paste0(opt\$out_prefix, "_matrix.csv"), quote=FALSE, row.names=TRUE)

    png(paste0(opt\$out_prefix, "_heatmap.png"), width=900, height=700)
    pheatmap(sim, display_numbers=TRUE, cluster_rows=FALSE, cluster_cols=FALSE)
    dev.off()

    # ---- Versions ----
    r_version <- paste0(R.version\$major, ".", R.version\$minor)
    writeLines(c(
        "\\"${task.process}\\":",
        paste0("    r-base: ", r_version),
        paste0("    r-lsa: ", packageVersion("lsa")),
        paste0("    r-pheatmap: ", packageVersion("pheatmap"))
    ), "versions.yml")
    """
}
