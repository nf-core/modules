#!/usr/bin/env Rscript

# Load necessary libraries
library(data.table)
library(stringr)
library(janitor)

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Check for Non-Empty, Non-Whitespace String
#'
#' This function checks if the input is non-NULL and contains more than
#' just whitespace.
#' It returns TRUE if the input is a non-empty, non-whitespace string,
#' and FALSE otherwise.
#'
#' @param input A variable to check.
#' @return A logical value: TRUE if the input is a valid, non-empty,
#' non-whitespace string; FALSE otherwise.
#' @examples
#' is_valid_string("Hello World") # Returns TRUE
#' is_valid_string("   ")         # Returns FALSE
#' is_valid_string(NULL)          # Returns FALSE
is_valid_string <- function(input) {
  !is.null(input) && nzchar(trimws(input))
}

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse
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

#' Turn “null” or empty strings into actual NULL
#'
#' @param x Input option
#'
#' @return NULL or x
#'
nullify <- function(x) {
  if (is.character(x) && (tolower(x) == "null" || x == "")) NULL else x
}

#' Parse tolerance value
#'
#' @param x Tolerance to check
#'
#' @return 1e-6 if x value is null and x
#' if x is > 0 and numeric else error
#'
#' @examples
#' parse_tolerance(NULL) # 1e-6
#' parse_tolerance("0.1") # 0.1
#' parse_tolerance(NA) # error
#' parse_tolerance("ABC") # error
#' parse_tolerance("-1.5") # error
parse_tolerance <- function(x) {
  if (is.null(x)) return(1e-6)
  out <- suppressWarnings(as.numeric(x))
  if (is.na(out)) stop("tolerance must be numeric")
  if (out < 0) stop("tolerance must be non-negative")
  out
}

#' Return the corresponding normalise column name
#'
#' @param x non normalised column name
#'
#' @return corresponding string if recognise or x
#'
#' @examples
#' convert_colnames(c(" chr", "Chromosome")
#' convert_colnames("genetic_map (cM).")
#' convert_colnames("Combined Rate (cM/Mb).")
convert_colnames <- function(x) {
  x <- janitor::make_clean_names(tolower(x))

  recode <- c(
    chr = "chr", `#chr` = "chr", chrom = "chr",
    chromosome = "chr", `_chr` = "chr",
    id = "id", snp = "id", marker = "id", rsid = "id",
    pos = "pos", position = "pos", bp = "pos",
    cm = "cm", genetic_map = "cm",
    genetic_map_cm = "cm", genetic_map_cm_ = "cm",
    rate = "rate", combined_rate = "rate", cm_mb = "rate",
    combined_rate_cm_mb_ = "rate", combined_rate_cm_mb = "rate"
  )

  ifelse(is.na(recode[x]), x, recode[x])
}

#' Main function to process the map file
#'
#' @param file_path Path to genetic map file to convert
#' @param chr Chromosome name
#' @param prefix Prefix name for the output files
#' @param tolerance Difference
process_map_file <- function(
  file_path, chr = NULL, prefix = "output",
  tolerance = NULL
) {
  # Read the map file into a data.table
  options(warn = 2) # all warnings will be set to error
  map_df <- data.table::fread(
    file_path,
    sep = "auto",
    header = "auto",
    showProgress = FALSE
  )

  if (nrow(map_df) == 0) stop("Input map is empty")

  no_header <- all(grepl(
    "^V[0-9]+\\\\z", colnames(map_df),
    ignore.case = TRUE, perl = TRUE
  ))
  if (no_header) {
    if (dim(map_df)[2] == 3) {
      message("Ambiguous no-header input, inferring to be chr, pos, cm")
      colnames(map_df) <- c("chr", "pos", "cm")
    } else if (dim(map_df)[2] == 4) {
      message("Ambiguous no-header input, inferring to be chr, id, cm, pos")
      colnames(map_df) <- c("chr", "id", "cm", "pos")
    } else {
      stop(
        "Cannot auto-detect column names for file with ",
        dim(map_df)[2],
        " columns without header"
      )
    }
  } else {
    colnames(map_df) <- convert_colnames(colnames(map_df))
  }

  # Initialize columns missing columns
  if (!"chr" %in% colnames(map_df)) {
    if (is.null(chr)) {
      stop("Chromosome column missing and chr not present in meta")
    }
    map_df[["chr"]] <- as.character(chr)
  } else {
    map_df[["chr"]] <- as.character(map_df[["chr"]])
  }

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
  map_df[, cm := cm - cm[1]]

  if (!isTRUE(all.equal(map_df[["cm"]][1], 0))) {
    stop("First cM value must be 0 after normalization")
  }

  # Compute forward rate for previous row (interval prev -> current)
  delta_bp <- c(diff(map_df[["pos"]]))
  delta_cm <- c(diff(map_df[["cm"]]))
  rate <- c(delta_cm / delta_bp * 1e6, 0)

  if (!"rate" %in% colnames(map_df)) {
    map_df[["rate"]] <- rate
  } else {
    map_df[["diff"]] <- abs(map_df[["rate"]] - rate)
    if (any(map_df[["diff"]] > tolerance)) {
      print(map_df[map_df[["diff"]] > tolerance, ])
      stop("cm[n] must equal cm[n-1] + ( (pos[n] - pos[n-1]) / 1e6 * rate[n-1])")
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


################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

# I've defined these in a single array like this so that we could go back to an
# optparse-driven method in future with module bin/ directories, rather than
# the template

# Set defaults and classes

opt <- list(
  output_prefix = "${prefix}",
  map_file = "${map_file}",
  chr = "${meta.chr}",
  tolerance = NULL
)

opt_types <- lapply(opt, class)

# Apply parameter overrides
args_opt <- parse_args("${args}")
for (ao in names(args_opt)) {
  if (! ao %in% names(opt)) {
    stop(paste("Invalid option:", ao))
  } else {
    # Preserve classes from defaults where possible
    if (! is.null(opt[[ao]])) {
      args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
    }
    opt[[ao]] <- args_opt[[ao]]
  }
}

keys <- c("tolerance", "chr")
opt[keys] <- lapply(opt[keys], nullify)

for (file_input in c("map_file")) {
  if (! is_valid_string(opt[[file_input]])) {
    stop(paste("Please provide", file_input), call. = FALSE)
  }

  if (! file.exists(opt[[file_input]])) {
    stop(paste0(
      "Value of ", file_input, ": ",
      opt[[file_input]], " is not a valid file"
    ))
  }
}

process_map_file(
  file_path = opt[["map_file"]],
  chr =  opt[["chr"]],
  prefix =  opt[["output_prefix"]],
  tolerance =  parse_tolerance(opt[["tolerance"]])
)

version_rbase <- paste(R.version[["major"]], R.version[["minor"]], sep = ".")
version_datatable <- packageVersion("data.table")
version_janitor <- packageVersion("janitor")

writeLines(c(
  '"${task.process}":',
  paste("    r-base:", version_rbase),
  paste("    r-data.table:", version_datatable),
  paste("    r-janitor:", version_janitor)
), "versions.yml")
