#!/usr/bin/env Rscript

he_overall_file <- "$he_overall_file"
he_within_file <- "$he_within_file"
he_across_file <- "$he_across_file"

parse_he_file <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(paste("File does not exist:", file_path))
    return(NULL)
  }

  lines <- readLines(file_path)
  her_all_line <- grep("^Her_All", lines, value = TRUE)

  if (length(her_all_line) == 0) {
    warning(paste("Her_All line not found in file:", file_path))
    return(NULL)
  }

  parts <- strsplit(her_all_line, "\\\\s+")[[1]]

  if (length(parts) < 3) {
    warning(paste("Invalid Her_All line format in file:", file_path))
    return(NULL)
  }

  data.frame(
    file = basename(file_path),
    heritability = as.numeric(parts[2]),
    se = as.numeric(parts[3]),
    stringsAsFactors = FALSE
  )
}

overall_result <- parse_he_file(he_overall_file)
within_result <- parse_he_file(he_within_file)
across_result <- parse_he_file(he_across_file)

h2_overall <- NA_real_
h2_overall_se <- NA_real_
h2_same <- NA_real_
h2_same_se <- NA_real_
h2_diff <- NA_real_
h2_diff_se <- NA_real_
T2_statistic <- NA_real_

statistical_test_results <- list(
  pvalue = NA_real_,
  mean_T2samp = NA_real_,
  sd_T2samp = NA_real_
)

if (!is.null(overall_result)) {
  h2_overall <- overall_result\$heritability
  h2_overall_se <- overall_result\$se
}

if (!is.null(within_result)) {
  h2_same <- within_result\$heritability
  h2_same_se <- within_result\$se
}

if (!is.null(across_result)) {
  h2_diff <- across_result\$heritability
  h2_diff_se <- across_result\$se
}

if (!is.na(h2_same) && !is.na(h2_diff)) {
  T2_statistic <- h2_same - h2_diff
}

if (!is.na(h2_same) && !is.na(h2_same_se) && !is.na(h2_diff) && !is.na(h2_diff_se)) {
  mean_t2 <- h2_same - h2_diff
  sd_t2 <- sqrt(h2_same_se^2 + h2_diff_se^2)

  statistical_test_results\$pvalue <- pnorm(0, mean = mean_t2, sd = sd_t2)
  statistical_test_results\$mean_T2samp <- mean_t2
  statistical_test_results\$sd_T2samp <- sd_t2
}

output_lines <- c(
  "LDAK Genotype Error Analysis Results (T2 Statistic)",
  "=====================================================",
  "",
  "Input Files:",
  paste("  Overall HE file:", basename(he_overall_file)),
  paste("  Within-batch HE file:", basename(he_within_file)),
  paste("  Across-batch HE file:", basename(he_across_file)),
  "",
  "Overall Results:",
  paste("  Heritability:", round(h2_overall, 6)),
  paste("  SE:", round(h2_overall_se, 6)),
  "",
  "Same-Batch Results (h2Same):",
  paste("  Heritability:", round(h2_same, 6)),
  paste("  SE:", round(h2_same_se, 6)),
  "",
  "Different-Batch Results (h2Diff):",
  paste("  Heritability:", round(h2_diff, 6)),
  paste("  SE:", round(h2_diff_se, 6)),
  "",
  "Genotype Error Analysis:",
  paste("  T2 Statistic (h2Same - h2Diff):", round(T2_statistic, 6)),
  "",
  "Interpretation:",
  "  T2 > 0: Same-batch samples are more similar than expected",
  "           (indicates potential genotyping errors or batch effects)",
  "  T2 ~= 0: No evidence of batch-related errors",
  "  T2 < 0: Unexpected pattern (may indicate over-correction)",
  "",
  "Statistical Test Results:",
  paste("  P-value (H0: T2 <= 0):", round(statistical_test_results\$pvalue, 6)),
  paste("  Mean T2 (analytical):", round(statistical_test_results\$mean_T2samp, 6)),
  paste("  SD T2 (analytical):", round(statistical_test_results\$sd_T2samp, 6)),
  "",
  "Significance:",
  ifelse(
    !is.na(statistical_test_results\$pvalue),
    ifelse(
      statistical_test_results\$pvalue < 0.05,
      "  SIGNIFICANT: Genotype errors detected (p < 0.05)",
      "  NOT SIGNIFICANT: No strong evidence of genotype errors (p >= 0.05)"
    ),
    "  Cannot determine (missing data)"
  )
)

writeLines(output_lines, "${meta.id}.txt")

writeLines(
  c(
    "\"${task.process}\":",
    paste("    r-base:", paste(R.version[['major']], R.version[['minor']], sep = "."))
  ),
  "versions.yml"
)
