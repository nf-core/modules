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

gwas <- fread("${sumstats}", data.table = FALSE)
hdl_ext_args <- parse_hdl_ext_args()

call_args <- c(
  list(
    gwas.df = gwas,
    LD.path = "${hdl_ref_panel_dir}",
    output.file = "${meta.id}.hdl.log"
  ),
  hdl_ext_args
)

result <- do.call(HDL::HDL.h2, call_args)

out <- data.frame(
  trait = "${meta.id}",
  h2 = as.numeric(result[["h2"]]),
  se = as.numeric(result[["h2.se"]]),
  p = as.numeric(result[["P"]]),
  eigen_use = as.character(result[["eigen.use"]]),
  stringsAsFactors = FALSE
)

fwrite(out, "${meta.id}.h2.tsv", sep = "\t", quote = FALSE, na = "NA")
