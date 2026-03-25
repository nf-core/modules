#!/usr/bin/env Rscript

# Load necessary libraries
library(data.table)
library(stringr)

# Function to parse external arguments
parse_args <- function(x) {
  args_list <- unlist(strsplit(x, " ?--")[[1]])[-1]
  args_vals <- lapply(
    args_list,
    function(x) scan(text = x, what = "character", quiet = TRUE)
  )

  # Ensure the option vectors are length 2 (key/ value) to catch empty ones
  args_vals <- lapply(args_vals, function(z) {
    length(z) <- 2
    z
  })

  parsed_args <- structure(
    lapply(args_vals, function(x) x[2]),
    names = lapply(args_vals, function(x) x[1])
  )
  parsed_args[! is.na(parsed_args)]
}

case_when_base <- function(x) {
  if (x %in% c("chr", "#chr", "chrom", "chromosome")) {
    "chr"
  } else if (x %in% c("id", "snp", "marker", "rsid")) {
    "id"
  } else if (x %in% c("pos", "position", "bp")) {
    "pos"
  } else if (x %in% c("cm", "genetic_map", "genetic_map_cm_")) {
    "cm"
  } else if (x %in% c(
    "rate", "combined_rate", "cm_cb", "cm_mb",
    "combined_rate_cm_mb_"
  )) {
    "rate"
  } else {
    x
  }
}

normalize_names <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9#]+", "_", x)

  sapply(x, case_when_base)
}

# Main function to process the map file
process_map_file <- function(
  file_path, chr = NULL, prefix = "output",
  tolerance = NA
) {
  if (chr == "null") {
    chr <- NULL
  }
  # Read the map file into a data.table
  options(warn = 2) # all warnings will be set to error
  map_df <- data.table::fread(
    file_path,
    sep = "auto",
    header = "auto",
    showProgress = FALSE
  )
  no_header <- all(stringr::str_detect(colnames(map_df), "^V[0-9]+\\\\z"))
  if (no_header) {
    if (dim(map_df)[2] == 3) {
      colnames(map_df) <- c("chr", "pos", "cm")
    } else if (dim(map_df)[2] == 4) {
      colnames(map_df) <- c("chr", "id", "cm", "pos")
    } else {
      stop(
        "Cannot auto-detect column names for file with ",
        dim(map_df)[2],
        " columns without header"
      )
    }
  } else {
    colnames(map_df) <- normalize_names(colnames(map_df))
  }

  # Initialize columns missing columns
  if (!"chr" %in% colnames(map_df)) {
    # Ensure chromosome column is present
    print(chr)
    if (is.null(chr)) {
      stop("Chromosome column missing and chr not present in meta")
    }
    map_df[["chr"]] <- as.character(chr)
  } else {
    map_df[["chr"]] <- as.character(map_df[["chr"]])
  }

  print(map_df)

  if (length(unique(map_df[["chr"]])) > 1) {
    stop("More than one chromosome present in file")
  }

  if (!"id" %in% colnames(map_df)) {
    map_df[["id"]] <- "."
  }

  # Ensure necessary columns are present
  if (!all(c("pos", "cm") %in% colnames(map_df))) {
    stop("Position and cM missing")
  }

  if (!is.numeric(map_df[["cm"]])) {
    stop("cM column should be numeric")
  }

  if (!is.numeric(map_df[["pos"]])) {
    stop("pos column should be numeric")
  }

  missing_row <- any(is.na(map_df[["pos"]]) | is.na(map_df[["cm"]]))

  if (missing_row) {
    stop("Position or cM missing")
  }

  if (!is.null(chr) && any(map_df[["chr"]] != chr)) {
    stop("Mismatch between chr given and the chr present in file")
  }

  # Order position to ensure all successive rows have increasing position
  map_df <- map_df[order(map_df[["pos"]]), ]
  if (any(duplicated(map_df[["pos"]]))) {
    print(map_df[duplicated(map_df[["pos"]]), ])
    stop("pos column shouldn't have any duplicate row")
  }

  # Normalize cM (needed by stitch)
  map_df[["cm"]] <- map_df[["cm"]] - map_df[["cm"]][1]

  if (map_df[["cm"]][1] != 0) {
    stop("cm[0] needs to be 0 for STITCH software")
  }

  # Compute forward rate for previous row (interval prev -> current)
  delta_bp <- c(diff(map_df[["pos"]]))
  delta_cm <- c(diff(map_df[["cm"]]))
  rate <- c(delta_cm / delta_bp * 1e6, 0)

  if (!"rate" %in% colnames(map_df)) {
    map_df[["rate"]] <- rate
  } else {
    map_df[["diff"]] <- abs(map_df[["rate"]] - rate)

    tolerance <- ifelse(
      is.na(tolerance) || is.null(tolerance),
      10e-6,
      as.numeric(tolerance)
    )
    if (any(map_df[["diff"]] > tolerance)) {
      print(map_df[map_df[["diff"]] > tolerance, ])
      stop("cm[n] must equal cm[n-1] + ( (pos[n] - pos[n-1]) * rate[n-1])")
    }
  }

  if (!is.numeric(map_df[["rate"]])) {
    stop("rate column should be numeric")
  }

  map_df[["rate"]] <- signif(map_df[["rate"]], digits = 8)
  map_df[["cm"]] <- signif(map_df[["cm"]], digits = 8)

  # Process the data
  glimpse_file <- paste0(prefix, ".glimpse.map")
  minimac_file <- paste0(prefix, ".minimac.map")
  plink_file <- paste0(prefix, ".plink.map")
  stitch_file <- paste0(prefix, ".stitch.map")

  # Write headers
  writeLines("pos\tchr\tcM", glimpse_file)
  writeLines("#chr\tposition\tGenetic_Map(cM)", minimac_file)
  writeLines("position COMBINED_rate.cM.Mb. Genetic_Map.cM.", stitch_file)

  # Write data to files
  fwrite(
    map_df[, c("pos", "chr", "cm")], glimpse_file,
    sep = "\t", append = TRUE, col.names = FALSE
  )
  fwrite(
    map_df[, c("chr", "pos", "cm")], minimac_file,
    sep = "\t", append = TRUE, col.names = FALSE
  )
  fwrite(
    map_df[, c("chr", "id", "cm", "pos")], plink_file,
    sep = " ", append = TRUE, col.names = FALSE
  )
  fwrite(
    map_df[, c("pos", "rate", "cm")], stitch_file,
    sep = " ", append = TRUE, col.names = FALSE
  )
}

ext_args <- parse_args("${args}")

process_map_file(
  "${map_file}",
  chr = "${meta.chr}",
  prefix = "${prefix}",
  tolerance = ext_args[["tolerance"]]
)

version_rbase <- paste(R.version[["major"]], R.version[["minor"]], sep = ".")
version_datatable <- packageVersion("data.table")
version_stringr <- packageVersion("stringr")

writeLines(c(
  '"${task.process}":',
  paste("    r-base:", version_rbase),
  paste("    r-data.table:", version_datatable),
  paste("    r-stringr:", version_stringr)
), "versions.yml")
