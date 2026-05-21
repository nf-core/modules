#!/usr/bin/env Rscript

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Check for Non-Empty, Non-Whitespace String
is_valid_string <- function(input) {
    !is.null(input) && nzchar(trimws(input))
}

#' Parse long-form options like --opt1 val1 --opt2 val2
parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})
    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

#' Flexibly read CSV / TSV / featureCounts-style tables
read_delim_flexible <- function(file, header = TRUE, sep = NULL, comment.char = "", check.names = TRUE){
    ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1))
    if (ext == "gz") {
        # peek at the inner extension
        inner <- tolower(tail(strsplit(sub("\\\\.gz\$", "", basename(file)), split = "\\\\.")[[1]], 1))
    } else {
        inner <- ext
    }
    if (is.null(sep)) {
        sep <- if (inner == "csv") "," else "\t"
    }
    read.table(
        file,
        sep = sep,
        header = header,
        comment.char = comment.char,
        stringsAsFactors = FALSE,
        check.names = check.names
    )
}

################################################
################################################
## Parse parameters from Nextflow             ##
################################################
################################################

opt <- list(
    output_prefix         = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    count_file            = '$counts',
    sample_file           = '$samplesheet',
    flattened_gtf         = '$flattened_gtf',
    flattened_bed         = '$flattened_bed',
    contrast_variable     = '$contrast_variable',
    reference_level       = '$reference',
    target_level          = '$target',
    sample_id_col         = "run",
    strategy_col          = "strategy",
    replicate_col         = "replicate",
    sample_name_regex     = NULL,
    modules               = "DOU,DTE",
    min_count             = as.integer(1),
    stringent             = "TRUE",
    dispersion_modeling   = "auto",
    nullweight            = as.numeric(500),
    contrasts_method      = "revpairwise",
    cores                 = as.integer('$task.cpus')
)
opt_types <- lapply(opt, class)

args_opt <- parse_args('$task.ext.args')
for (ao in names(args_opt)) {
    if (!ao %in% names(opt)) stop(paste("Invalid option:", ao))
    if (!is.null(opt[[ao]])) {
        args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
    }
    opt[[ao]] <- args_opt[[ao]]
}

required_opts <- c("contrast_variable", "reference_level", "target_level", "output_prefix")
missing <- required_opts[!unlist(lapply(opt[required_opts], is_valid_string)) | !required_opts %in% names(opt)]
if (length(missing) > 0) {
    stop(paste("Missing required options:", paste(missing, collapse = ", ")))
}

for (file_input in c("count_file", "sample_file", "flattened_gtf", "flattened_bed")){
    if (!is_valid_string(opt[[file_input]])) stop(paste("Please provide", file_input))
    if (!file.exists(opt[[file_input]])) stop(paste0("Value of ", file_input, ": ", opt[[file_input]], " is not a valid file"))
}

modules <- trimws(strsplit(opt\$modules, ",")[[1]])
stringent_val <- switch(toupper(opt\$stringent),
    "TRUE" = TRUE,
    "FALSE" = FALSE,
    "NULL" = NULL,
    stop("`stringent` must be one of TRUE, FALSE, NULL")
)

################################################
################################################
## Load libraries                             ##
################################################
################################################

suppressPackageStartupMessages({
    library(DOTSeq)
    library(SummarizedExperiment)
})

################################################
################################################
## Read inputs                                ##
################################################
################################################

cnt <- read_delim_flexible(opt\$count_file, header = TRUE, comment.char = "#")

# Allow optional renaming of long count-table column names via regex
if (!is.null(opt\$sample_name_regex) && nzchar(opt\$sample_name_regex)) {
    names(cnt) <- gsub(opt\$sample_name_regex, "\\\\1", names(cnt))
}

cond <- read_delim_flexible(opt\$sample_file, header = TRUE)

# Normalise condition column names and rename user-chosen columns into the
# names that DOTSeq expects ("run", "strategy", "replicate", "condition")
names(cond) <- tolower(trimws(names(cond)))
user_to_required <- list(
    run = tolower(opt\$sample_id_col),
    strategy = tolower(opt\$strategy_col),
    replicate = tolower(opt\$replicate_col),
    condition = tolower(opt\$contrast_variable)
)

for (req in names(user_to_required)) {
    src <- user_to_required[[req]]
    if (!src %in% names(cond)) {
        stop(paste0("Sample sheet column '", src, "' (mapped to '", req, "') not found. Have: ",
                    paste(names(cond), collapse = ", ")))
    }
    if (src != req) {
        # Drop any pre-existing column with the required name to avoid collision
        cond[[req]] <- cond[[src]]
        if (src %in% names(cond) && src != req) cond[[src]] <- NULL
    }
}

# Filter samplesheet to only the levels involved in the contrast
cond <- cond[cond\$condition %in% c(opt\$reference_level, opt\$target_level), , drop = FALSE]
if (nrow(cond) == 0) {
    stop("No samples remain after filtering condition column to reference/target levels.")
}

# Set baseline as the reference level so DTE coefficient is target_vs_reference
cond\$condition <- factor(cond\$condition, levels = c(opt\$reference_level, opt\$target_level))

################################################
################################################
## Build DOTSeqDataSets                       ##
################################################
################################################

d <- DOTSeqDataSetsFromFeatureCounts(
    count_table = cnt,
    condition_table = cond,
    flattened_gtf = opt\$flattened_gtf,
    flattened_bed = opt\$flattened_bed,
    min_count = opt\$min_count,
    stringent = stringent_val,
    baseline = opt\$reference_level,
    verbose = FALSE
)

################################################
################################################
## Run DOTSeq                                 ##
################################################
################################################

d <- DOTSeq(
    datasets = d,
    modules = modules,
    target = opt\$target_level,
    baseline = opt\$reference_level,
    min_count = opt\$min_count,
    stringent = stringent_val,
    dispersion_modeling = opt\$dispersion_modeling,
    nullweight = opt\$nullweight,
    contrasts_method = opt\$contrasts_method,
    parallel = list(n = opt\$cores, autopar = TRUE),
    verbose = FALSE
)

################################################
################################################
## Extract results                            ##
################################################
################################################

write_results <- function(df, suffix) {
    if (is.null(df) || (is.data.frame(df) && nrow(df) == 0)) return(invisible(NULL))
    out_df <- as.data.frame(df)
    write.table(
        out_df,
        file = paste(opt\$output_prefix, suffix, sep = "."),
        col.names = TRUE,
        row.names = FALSE,
        sep = "\t",
        quote = FALSE
    )
}

interaction_results <- tryCatch(getContrasts(d, type = "interaction"), error = function(e) NULL)
if (!is.null(interaction_results)) {
    write_results(interaction_results\$DOU, "dou.interaction.dotseq.results.tsv")
    write_results(interaction_results\$DTE, "dte.interaction.dotseq.results.tsv")
}

strategy_results <- tryCatch(getContrasts(d, type = "strategy"), error = function(e) NULL)
if (!is.null(strategy_results)) {
    write_results(strategy_results\$DOU, "dou.strategy.dotseq.results.tsv")
    write_results(strategy_results\$DTE, "dte.strategy.dotseq.results.tsv")
}

# Serialise the full DOTSeqDataSets object for downstream use
saveRDS(d, file = paste(opt\$output_prefix, "DOTSeqDataSets.rds", sep = "."))

################################################
################################################
## R session info                             ##
################################################
################################################

sink(paste(opt\$output_prefix, "R_sessionInfo.log", sep = "."))
print(sessionInfo())
sink()

################################################
################################################
## Versions                                   ##
################################################
################################################

dotseq.version <- as.character(packageVersion("DOTSeq"))

writeLines(
    c(
        '"${task.process}":',
        paste("    bioconductor-dotseq:", dotseq.version)
    ),
    "versions.yml"
)
