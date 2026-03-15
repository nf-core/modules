#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

parse_hdl_ext_args <- function(raw_ext_args = "${task.ext.args ?: ''}") {
  tryCatch(
    eval(parse(text = sprintf("list(%s)", raw_ext_args))),
    error = function(e) stop(sprintf("Failed to parse task.ext.args as R named arguments: %s", conditionMessage(e)))
  )
}

if (!requireNamespace("HDL", quietly = TRUE)) {
  stop("R package 'HDL' is not installed. Install HDL in the runtime container/environment.")
}

gwas1 <- fread("${sumstats1}", data.table = FALSE)
gwas2 <- fread("${sumstats2}", data.table = FALSE)
hdl_ext_args <- parse_hdl_ext_args()

call_args <- c(
  list(
    gwas1.df = gwas1,
    gwas2.df = gwas2,
    LD.path = "${hdl_ref_panel_dir}",
    output.file = "${meta.id}.${meta2.id}.hdl.log"
  ),
  hdl_ext_args
)

result <- do.call(HDL::HDL.rg, call_args)

h1 <- NA_real_
h2 <- NA_real_
gcov <- NA_real_
if (!is.null(result[["estimates.df"]])) {
  est <- result[["estimates.df"]]
  if ("Heritability_1" %in% rownames(est)) {
    h1 <- as.numeric(est["Heritability_1", "Estimate"])
  }
  if ("Heritability_2" %in% rownames(est)) {
    h2 <- as.numeric(est["Heritability_2", "Estimate"])
  }
  if ("Genetic_Covariance" %in% rownames(est)) {
    gcov <- as.numeric(est["Genetic_Covariance", "Estimate"])
  }
}

out <- data.frame(
  trait1 = "${meta.id}",
  trait2 = "${meta2.id}",
  rg = as.numeric(result[["rg"]]),
  se = as.numeric(result[["rg.se"]]),
  p = as.numeric(result[["P"]]),
  h2_trait1 = h1,
  h2_trait2 = h2,
  covariance = gcov,
  eigen_use = as.character(result[["eigen.use"]]),
  stringsAsFactors = FALSE
)

fwrite(out, "${meta.id}.${meta2.id}.rg.tsv", sep = "\t", quote = FALSE, na = "NA")
