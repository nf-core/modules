#!/usr/bin/env Rscript

################################################
## Run CoGAPS (Bayesian non-negative matrix
## factorisation) on a features x samples matrix
################################################

## Options are passed via task.ext.args as '--name value' pairs, where name is
## a scalar CogapsParams parameter, e.g. '--nPatterns 7 --nIterations 50000 --seed 42'.
## Non-scalar parameters (e.g. fixedPatterns, explicitSets, subsetIndices,
## geneNames) cannot be set this way. See ?CogapsParams for the full list.

suppressPackageStartupMessages(library(CoGAPS))

parse_args <- function(x) {
    x <- trimws(x)
    if (x == "" || x == "null") return(list())
    parts <- trimws(strsplit(x, "--", fixed = TRUE)[[1]])
    parts <- parts[parts != ""]
    out <- list()
    for (p in parts) {
        kv <- strsplit(p, "[[:space:]=]+")[[1]]
        out[[kv[1]]] <- if (length(kv) > 1) paste(kv[-1], collapse = " ") else "TRUE"
    }
    out
}

read_matrix <- function(file) {
    ext <- tolower(tools::file_ext(sub("[.]gz", "", file)))
    if (ext == "rds") {
        obj <- readRDS(file)
        if (is.data.frame(obj)) obj <- as.matrix(obj)
        return(obj)
    }
    if (ext == "csv") {
        df <- read.csv(file, row.names = 1, check.names = FALSE)
    } else if (ext %in% c("tsv", "txt")) {
        df <- read.delim(file, row.names = 1, check.names = FALSE)
    } else {
        stop("Unsupported input format '", ext, "': use csv, tsv, txt (optionally gzipped) or rds")
    }
    as.matrix(df)
}

prefix <- if ('${task.ext.prefix}' == 'null') '${meta.id}' else '${task.ext.prefix}'
opt    <- parse_args('${task.ext.args}')

## nPatterns has no default in CoGAPS, so it must be given
if (is.null(opt[["nPatterns"]])) {
    stop("nPatterns is required: set it in ext.args, e.g. '--nPatterns 7'")
}

## Build CoGAPS parameters from ext.args, coercing each value to the slot's type
params <- CogapsParams(nPatterns = as.integer(opt[["nPatterns"]]))
for (name in setdiff(names(opt), "nPatterns")) {
    if (!name %in% slotNames(params)) {
        stop("Unknown CoGAPS parameter '", name, "'. Valid: ", paste(slotNames(params), collapse = ", "))
    }
    value <- opt[[name]]
    type <- class(slot(params, name))[1]
    coerced <- suppressWarnings(switch(type,
        numeric   = as.numeric(value),
        integer   = as.integer(value),
        logical   = as.logical(value),
        value))
    if (length(coerced) != 1 || is.na(coerced)) {
        stop("Invalid value '", value, "' for CoGAPS parameter '", name, "' (expected ", type, ")")
    }
    params <- setParam(params, name, coerced)
}

data <- read_matrix('${matrix}')

result <- CoGAPS(
    data,
    params   = params,
    nThreads = as.integer('${task.cpus}'),
    messages = TRUE
)

saveRDS(result, paste0(prefix, ".cogaps.rds"))

loadings <- getFeatureLoadings(result)
factors  <- getSampleFactors(result)
write.csv(data.frame(feature = rownames(loadings), loadings, check.names = FALSE),
          paste0(prefix, ".feature_loadings.csv"), row.names = FALSE)
write.csv(data.frame(sample = rownames(factors), factors, check.names = FALSE),
          paste0(prefix, ".sample_factors.csv"), row.names = FALSE)

## versions.yml
writeLines(c(
    '"${task.process}":',
    paste0("    cogaps: ", as.character(packageVersion("CoGAPS"))),
    paste0("    r-base: ", as.character(getRversion()))
), "versions.yml")
