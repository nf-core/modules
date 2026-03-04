#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(lsa)
  library(pheatmap)
})

# ---- Configurable options via task.ext.args ----
option_list <- list(
  make_option(c("-m", "--method"),
    type = "character", default = "cosine",
    help = "Similarity method: cosine, pearson, spearman"),
  make_option(c("-g", "--min_gene_mean"),
    type = "double", default = 0.0,
    help = "Minimum gene mean expression filter")
)

opt <- parse_args(
  OptionParser(option_list = option_list),
  args = strsplit("$args", "[[:space:]]+")[[1]]
)

# ---- Load table ----
df <- as.data.frame(
  read_csv("${expression_matrix}",
    col_types = cols(.default = col_guess()))
)

# ---- Set rownames from first column if character ----
first_col <- names(df)[1]
if (is.character(df[[first_col]]) || is.factor(df[[first_col]])) {
  rownames(df) <- df[[first_col]]
  df[[first_col]] <- NULL
}

# ---- Select numeric columns as the expression matrix ----
mat <- as.matrix(
  df[, sapply(df, is.numeric), drop = FALSE]
)

# ---- Filter low-expression genes ----
if (opt\$min_gene_mean > 0) {
  keep <- rowMeans(mat, na.rm = TRUE) >= opt\$min_gene_mean
  mat <- mat[keep, , drop = FALSE]
}

# ---- Compute similarity ----
if (opt\$method == "cosine") {
  sim <- lsa::cosine(mat)
} else {
  sim <- cor(mat,
    method = opt\$method, use = "pairwise.complete.obs")
}
colnames(sim) <- colnames(mat)
rownames(sim) <- colnames(mat)

# ---- Save matrix ----
write.csv(sim, "${prefix}_matrix.csv",
  quote = FALSE, row.names = TRUE)

# ---- Save heatmap ----
png("${prefix}_heatmap.png", width = 900, height = 700)
pheatmap(sim,
  display_numbers = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE)
dev.off()

# ---- Versions ----
writeLines(
  c(
    '"${task.process}":',
    paste("    r-base:",
      paste0(R.version\$major, ".", R.version\$minor)),
    paste("    r-lsa:",
      as.character(packageVersion("lsa"))),
    paste("    r-pheatmap:",
      as.character(packageVersion("pheatmap")))
  ),
  "versions.yml"
)
