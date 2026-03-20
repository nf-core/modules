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

read_first_line <- function(file_path) {
  is_url <- grepl("^https?://", file_path)
  con <- if (is_url) url(file_path, "rb") else file(file_path, "rb")

  if (endsWith(file_path, ".gz")) {
    con <- gzcon(con)
  }

  on.exit(close(con), add = TRUE)
  readLines(con, n = 1, warn = FALSE)
}

# Function to detect separator
detect_separator <- function(first_line) {
  counts <- c(
    tab = str_count(first_line, "\\t"),
    comma = str_count(first_line, ","),
    semicolon = str_count(first_line, ";"),
    space = str_count(first_line, " ")
  )

  separators <- c("\\t", ",", ";", " ")
  separators[which.max(counts)]
}

# Function to detect header
detect_header <- function(first_line) {
  header_keywords <- c(
    "pos", "cm", "snp",
    "position", "Genetic_Map",
    "rate", "COMBINED_rate"
  )
  if (any(sapply(
    header_keywords, function(keyword) {
      grepl(keyword, first_line, ignore.case = TRUE)
    }
  ))) {
    TRUE
  } else {
    FALSE
  }
}

# Function to detect column names
detect_column_names <- function(first_line, detected_sep, detected_header) {
  if (detected_header) {
    cols <- str_split(first_line, detected_sep)[[1]]
    cols <- tolower(cols)
    matching_cols <- list(
      "combined_rate" = "rate",
      "cm.cb" = "rate",
      "cm/mb" = "rate",
      "position" = "pos",
      "genetic_map" = "cm",
      "#chr" = "chr"
    )
    sapply(cols, function(col) {
      # Then check pattern matches
      for (pattern in names(matching_cols)) {
        if (grepl(pattern, col, fixed = TRUE)) {
          return(matching_cols[[pattern]])
        }
      }
      col
    })
  } else {
    num_cols <- length(str_split(first_line, detected_sep)[[1]])
    if (num_cols == 3) {
      c("chr", "pos", "cm")
    } else if (num_cols == 4) {
      c("chr", "id", "cm", "pos")
    } else {
      stop(
        "Error: Cannot auto-detect column names for file with",
        num_cols,
        " columns without header"
      )
    }
  }
}

# Main function to process the map file
process_map_file <- function(
  file_path, chr = NULL, prefix = "output",
  tolerance = NA
) {
  # Read the first line
  first_line <- read_first_line(file_path)
  cat(first_line, "\n")

  # Check if first line is empty
  if (str_trim(first_line) == "") {
    stop("Error: First line is empty")
  }

  # Auto-detect format
  detected_sep <- detect_separator(first_line)
  cat("Detected: SEP='", detected_sep, "'\n")
  detected_header <- detect_header(first_line)
  cat("Detected: HEADER=", detected_header, "\n")
  detected_cols <- detect_column_names(
    first_line, detected_sep, detected_header
  )
  cat("Detected: COLS=", detected_cols, "\n")

  # Read the map file into a data.table
  map_df <- fread(file_path, sep = detected_sep, header = detected_header)
  print(map_df)
  colnames(map_df) <- detected_cols

  # Initialize columns missing columns
  if (!"chr" %in% colnames(map_df)) {
    # Ensure chromosome column is present
    if (is.null(chr)) {
      stop("Error: Chromosome column missing and chr not present in meta")
    }
    map_df[["chr"]] <- as.character(chr)
  } else {
    map_df[["chr"]] <- as.character(map_df[["chr"]])
  }

  if (!"id" %in% colnames(map_df)) {
    map_df[["id"]] <- "."
  }

  # Ensure necessary columns are present
  if (!all(c("pos", "cm") %in% colnames(map_df))) {
    stop("Error: Position and cM missing")
  }

  if (!is.numeric(map_df[["cm"]])) {
    stop("Error: cM column should be numeric")
  }

  if (!is.numeric(map_df[["pos"]])) {
    stop("Error: pos column should be numeric")
  }

  missing_row <- any(is.na(map_df[["pos"]]) | is.na(map_df[["cm"]]))

  if (missing_row) {
    stop("Error: Position or cM missing")
  }

  if (!is.null(chr) && any(map_df[["chr"]] != chr)) {
    stop("Error: mismatch between chr given and the chr present in file")
  }

  # Compute forward rate for previous row (interval prev -> current)
  delta_bp <- c(diff(map_df[["pos"]]), NA)
  delta_cm <- c(diff(map_df[["cm"]]), NA)
  rate <- ifelse(
    delta_bp > 0 & !is.na(delta_bp),
    (delta_cm / delta_bp * 1e6), 0
  )

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
    stop("Error: pos column should be numeric")
  }

  # Process the data
  glimpse_file <- paste0(prefix, ".glimpse.map")
  minimac_file <- paste0(prefix, ".minimac.map")
  plink_file <- paste0(prefix, ".plink.map")
  stitch_file <- paste0(prefix, ".stitch.map")
  eagle_file <- paste0(prefix, ".eagle.map")

  # Write headers
  writeLines("pos\tchr\tcM", glimpse_file)
  writeLines("#chr\tposition\tGenetic_Map(cM)", minimac_file)
  writeLines("position COMBINED_rate.cM.Mb. Genetic_Map.cM.", stitch_file)
  writeLines("chr position COMBINED_rate.cM.Mb. Genetic_Map.cM.", eagle_file)

  # Write data to files
  con_glimpse <- file(glimpse_file, open = "a")
  con_minimac <- file(minimac_file, open = "a")
  con_plink <- file(plink_file, open = "a")
  con_stitch <- file(stitch_file, open = "a")
  con_eagle <- file(stitch_file, open = "a")

  # Write data to files
  writeLines(
    paste(map_df[["pos"]], map_df[["chr"]], map_df[["cm"]], sep = "\t"),
    con_glimpse
  )
  writeLines(
    paste(map_df[["chr"]], map_df[["pos"]], map_df[["cm"]], sep = "\t"),
    con_minimac
  )
  writeLines(
    paste(
      map_df[["chr"]], map_df[["id"]],
      map_df[["cm"]], map_df[["pos"]],
      sep = " "
    ),
    con_plink
  )
  writeLines(
    paste(map_df[["pos"]], map_df[["rate"]], map_df[["cm"]], sep = " "),
    con_stitch
  )
  writeLines(
    paste(
      map_df[["chr"]], map_df[["pos"]],
      map_df[["rate"]], map_df[["cm"]], sep = " "
    ),
    con_eagle
  )
  close(con_glimpse)
  close(con_minimac)
  close(con_plink)
  close(con_stitch)
  close(con_eagle)
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
