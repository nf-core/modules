#!/usr/bin/env Rscript

ldak_reml_file <- "$ldak_reml_file"
quarter_files <- c($quarter_reml_files_r)
output_file <- "$output_file"

parse_reml_file <- function(file_path) {
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

quarter_results <- list()
for (file in quarter_files) {
  result <- parse_reml_file(file)
  if (!is.null(result)) {
    result\$type <- "quarter"
    quarter_results[[length(quarter_results) + 1]] <- result
  }
}

ldak_results <- list()
ldak_result <- parse_reml_file(ldak_reml_file)
if (!is.null(ldak_result)) {
  ldak_result\$type <- "ldak"
  ldak_results[[1]] <- ldak_result
}

all_results <- do.call(rbind, c(quarter_results, ldak_results))
if (is.null(all_results) || nrow(all_results) == 0) {
  stop("No valid REML results found")
}

quarter_data <- all_results[all_results\$type == "quarter", , drop = FALSE]
ldak_data <- all_results[all_results\$type == "ldak", , drop = FALSE]

quarter_mean_h2 <- if (nrow(quarter_data) > 0) mean(quarter_data\$heritability, na.rm = TRUE) else NA_real_
quarter_mean_se <- if (nrow(quarter_data) > 0) mean(quarter_data\$se, na.rm = TRUE) else NA_real_
ldak_h2 <- if (nrow(ldak_data) > 0) ldak_data\$heritability[1] else NA_real_
ldak_se <- if (nrow(ldak_data) > 0) ldak_data\$se[1] else NA_real_

inflation_T1 <- NA_real_
inflation_factor <- NA_real_

if (!is.na(quarter_mean_h2) && !is.na(ldak_h2) && nrow(quarter_data) > 0) {
  n_quarters <- nrow(quarter_data)
  quarter_sum_h2 <- sum(quarter_data\$heritability, na.rm = TRUE)
  inflation_T1 <- (quarter_sum_h2 - ldak_h2) / (n_quarters - 1)
  inflation_factor <- quarter_mean_h2 / ldak_h2
}

statistical_test_results <- list(
  n_quarters = 0,
  pvalue = NA_real_,
  mean_T1samp = NA_real_,
  sd_T1samp = NA_real_
)

if (!is.na(ldak_h2) && !is.na(ldak_se) && nrow(quarter_data) > 0) {
  quarter_h2_values <- quarter_data\$heritability
  quarter_se_values <- quarter_data\$se

  valid_quarters <- !is.na(quarter_h2_values) & !is.na(quarter_se_values)
  quarter_h2_values <- quarter_h2_values[valid_quarters]
  quarter_se_values <- quarter_se_values[valid_quarters]

  if (length(quarter_h2_values) > 1) {
    mean_t1 <- (sum(quarter_h2_values) - ldak_h2) / (length(quarter_h2_values) - 1)
    sd_t1 <- sqrt(sum(quarter_se_values^2) + ldak_se^2) / (length(quarter_h2_values) - 1)

    statistical_test_results\$n_quarters <- length(quarter_h2_values)
    statistical_test_results\$pvalue <- pnorm(0, mean = mean_t1, sd = sd_t1)
    statistical_test_results\$mean_T1samp <- mean_t1
    statistical_test_results\$sd_T1samp <- sd_t1
  }
}

output_lines <- c(
  "LDAK Inflation Analysis Results",
  "================================",
  "",
  paste("Number of quarter files processed:", length(quarter_files)),
  paste("LDAK file processed:", basename(ldak_reml_file)),
  "",
  "Quarter Results:",
  paste("  Mean Heritability:", round(quarter_mean_h2, 6)),
  paste("  Mean SE:", round(quarter_mean_se, 6)),
  "",
  "LDAK Results:",
  paste("  Heritability:", round(ldak_h2, 6)),
  paste("  SE:", round(ldak_se, 6)),
  "",
  "Inflation Analysis:",
  paste("  T1 Statistic (Documentation Formula):", round(inflation_T1, 6)),
  paste("  Inflation Ratio (Legacy, Quarter/LDAK):", round(inflation_factor, 6)),
  "  Interpretation:",
  "    - T1 ~= 0: No inflation (good quality control)",
  "    - T1 > 0.05: Possible population structure or relatedness not fully captured",
  "    - T1 < -0.05: Possible over-correction",
  "",
  "Statistical Test Results:",
  paste("  Number of quarters used:", statistical_test_results\$n_quarters),
  paste("  P-value:", round(statistical_test_results\$pvalue, 6)),
  paste("  Mean T1 (analytical):", round(statistical_test_results\$mean_T1samp, 6)),
  paste("  SD T1 (analytical):", round(statistical_test_results\$sd_T1samp, 6)),
  "",
  "Individual Results:"
)

for (i in seq_len(nrow(all_results))) {
  row <- all_results[i, ]
  output_lines <- c(
    output_lines,
    paste(
      " ", row\$file, "(", row\$type, "):",
      "H2 =", round(row\$heritability, 6),
      "SE =", round(row\$se, 6)
    )
  )
}

writeLines(output_lines, output_file)

writeLines(
  c(
    "\"${task.process}\":",
    paste("    r-base:", paste(R.version[['major']], R.version[['minor']], sep = "."))
  ),
  "versions.yml"
)
