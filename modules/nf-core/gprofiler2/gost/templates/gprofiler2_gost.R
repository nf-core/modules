#!/usr/bin/env Rscript

# Written by Oskar Wacker (https://github.com/WackerO) in
# collaboration with Gisela Gabernet (https://github.com/ggabernet)
# Script template by Jonathan Manning (https://github.com/pinin4fjords)

# MIT License

# Copyright (c) QBiC

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x) {
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z) { length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = F) {

    ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1))

    if (ext == "tsv" || ext == "txt") {
        separator <- "\\t"
    } else if (ext == "csv") {
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }

    read.delim(
        file,
        sep = separator,
        header = header,
        row.names = row.names,
        check.names = check.names
    )
}

#' Round numeric dataframe columns to fixed decimal places by applying
#' formatting and converting back to numerics
#'
#' @param dataframe A data frame
#' @param columns Which columns to round (assumes all of them by default)
#' @param digits How many decimal places to round to? If -1, will return the unchanged input df
#'
#' @return output Data frame
round_dataframe_columns <- function(df, columns = NULL, digits = -1) {
    if (digits == -1) {
        return(df)                              # if -1, return df without rounding
    }

    df <- data.frame(df, check.names = FALSE)   # make data.frame from vector as otherwise, the format will get messed up
    if (is.null(columns)) {
        columns <- colnames(df)
    }
    df[,columns] <- round(
        data.frame(df[, columns], check.names = FALSE),
        digits = digits
    )

    # Convert columns back to numeric

    for (c in columns) {
        df[[c]][grep("^ *NA\$", df[[c]])] <- NA
        df[[c]] <- as.numeric(df[[c]])
    }
    df
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
    output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    de_file = '$de_file',
    de_id_column = 'gene_id',
    organism = NULL,
    significant = T,
    measure_underrepresentation = F,
    correction_method = 'gSCS',
    sources = NULL,
    evcodes = F,
    pval_threshold = '0.05',
    gmt_file = NULL,
    gost_id = NULL,
    background_file = '$background_file',
    domain_scope = 'annotated',
    min_diff = 1,
    round_digits = -1,
    enrich_colors = '#132B43,#56B1F7'
)

opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)) {
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

# Check if required parameters have been provided

required_opts <- c('output_prefix', 'organism', 'sources')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]

if (length(missing) > 0) {
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid

for (file_input in c('de_file')) {
    if (is.null(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }

    if (! file.exists(opt[[file_input]])) {
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(gprofiler2)
library(ggplot2)

################################################
################################################
## READ IN DIFFERENTIAL GENES FILE            ##
################################################
################################################

de.genes <-
    read_delim_flexible(
        file = opt\$de_file
    )
if (nrow(de.genes) == 0) {
    opt <- options(show.error.messages = FALSE) # Do not throw error so that rest of pipeline continues
    print("No differential features found, pathway enrichment analysis with gprofiler2 will be skipped.")
    stop()
}
query <- de.genes[[opt\$de_id_column]]

################################################
################################################
# Run gprofiler processes and generate outputs #
################################################
################################################

set.seed(1) # This will ensure that reruns have the same plot colors

output_prefix <- opt\$output_prefix
dir.create(output_prefix)

# Create empty output table in case no enriched pathways are found
file.create(paste(output_prefix, 'gprofiler2', 'all_enriched_pathways_tab', 'tsv', sep = '.'))

sources <-  strsplit(opt\$sources, split = ",")[[1]]
enrich_colors <- strsplit(opt\$enrich_colors, split = ",")[[1]]

if (!is.null(opt\$gost_id)) {

#   First check if a gost_id was provided
    gost_id <- opt\$gost_id
} else if (!is.null(opt\$gmt_file)){

    # Next check if custom GMT file was provided; extract only requested datasources (gprofiler will NOT filter automatically!)
    gmt <- Filter(function(line) any(startsWith(line, sources)), readLines(opt\$gmt))
    gmt_path <- paste0(strsplit(basename(opt\$gmt_file), split = "\\\\.")[[1]][[1]], ".", paste(sources, collapse="_"), "_filtered.gmt")
    writeLines(gmt, gmt_path)
    gost_id <- upload_GMT_file(gmt_path)
    
    # Add gost ID to output GMT name so that it can be reused in future runs
    file.rename(gmt_path, paste0(strsplit(basename(opt\$gmt_file), split = "\\\\.")[[1]][[1]], ".", paste(sources, collapse="_"), "_gostID_", gost_id, "_filtered.gmt"))
} else {

    # Otherwise, get the GMT file from gprofiler and save both the full file as well as the filtered one to metadata
    # TODO: Should this GMT also be filtered?
    gmt_url <- paste0("https://biit.cs.ut.ee/gprofiler//static/gprofiler_full_", opt\$organism, ".ENSG.gmt")
    tryCatch(
        {
            wget_command <- paste0("wget ", gmt_url)
            sys_return <- system(wget_command)
            if (sys_return != 0 && !("saved" %in% sys_return)) {
                print("Failed to fetch the GMT file from gprofiler with this URL:")
                print(gmt_url)
                print("For reproducibility reasons, try to download the GMT file manually by visiting https://biit.cs.ut.ee/gprofiler/gost, then selecting the correct organism and, in datasources, clicking 'combined ENSG.gmt'.")
            } else {
                gmt_path <- paste0("gprofiler_full_", opt\$organism, ".ENSG.gmt")
                gmt <- Filter(function(line) any(startsWith(line, sources)), readLines(gmt_path))
                gmt_path <- paste0("gprofiler_full_", opt\$organism, ".", paste(sources, collapse="_"), ".ENSG_filtered.gmt")
                writeLines(gmt, gmt_path)
            }
        },
        error=function(gost_error) {
            print("Failed to fetch the GMT file from gprofiler with this URL:")
            print(gmt_url)
            print("Got error:")
            print(gost_error)
            print("For reproducibility reasons, please try to download the GMT file manually by visiting https://biit.cs.ut.ee/gprofiler/gost, then selecting the correct organism and, in datasources, clicking 'combined ENSG.gmt'. Then provide it to the pipeline with the parameter `--gmt_file`")
        }
    )
    gost_id <- opt\$organism
}

capture.output(opt\$background_file, file="/home-link/iivow01/git/modules/error_gpro/ext_null")
capture.output("is.null(optbackground_file)", file="/home-link/iivow01/git/modules/error_gpro/ext_null", append=T)
capture.output((""==opt\$background_file), file="/home-link/iivow01/git/modules/error_gpro/ext_null", append=T)

# If custom background_file was provided, read it
if (opt\$background_file != "") {
    ext <- basename(opt\$background_file)
    capture.output(ext, file="/home-link/iivow01/git/modules/error_gpro/ext")
    ext <- strsplit(basename(opt\$background_file), split = "\\\\.")
    capture.output(ext, file="/home-link/iivow01/git/modules/error_gpro/ext", append=T)
#    ext <- tolower(tail(strsplit(basename(opt\$background_file), split = "\\\\.")[[1]], 1))
    if (ext == "txt") {
        background <- readLines(opt\$background_file)
    }

    # Check if tsv/csv or if second line contains tab or comma (in txt files); if so, the file is a table
    if (ext %in% c("csv", "tsv") || grepl("\\t|,", background[2]))
        intensities_table <- read_delim_flexible(
        file = opt\$background_file,
        row.names = 1
    )

    # Keep only numeric columns
    nums <- unlist(lapply(intensities_table, is.numeric), use.names = FALSE)
    intensities_table <- intensities_table[, nums]

    # Keep only rownames which have abundance
    background <- rownames(subset(intensities_table, rowSums(intensities_table, na.rm = TRUE)>0))
} else {
    background <- NULL
}

gost_results <- gost(
    query=query,
    organism=gost_id,
    significant=opt\$significant,
    measure_underrepresentation=opt\$measure_underrepresentation,
    correction_method=opt\$correction_method,
    sources=sources,
    evcodes=opt\$evcodes,
    user_threshold=opt\$pval_threshold,
    custom_bg=background,
    domain_scope=opt\$domain_scope
)

if (!is.null(gost_results)) {

    # Create interactive plot and save to HTML
    interactive_plot <- gostplot(gost_results, capped=T, interactive=T)

    # Save interactive plot as HTML
    htmlwidgets::saveWidget(
                widget = interactive_plot,
                file = paste(output_prefix, 'gprofiler2', 'gostplot', 'html', sep = '.')
                )

    # Create a static plot and save to PNG
    static_plot <- gostplot(gost_results, capped=T, interactive=F)
    ggsave(plot = static_plot, filename = paste(output_prefix, 'gprofiler2', 'gostplot', 'png', sep = '.'), width = 10, height = 7)

    # Subset gost results to those pathways with a min. number of differential features
    gost_results\$result <- gost_results\$result[which(gost_results\$result\$intersection_size>=opt\$min_diff),]

    # annotate query size (number of differential features in contrast)
    gost_results\$result\$original_query_size <- rep(length(as.character(de.genes\$Ensembl_ID)), nrow(gost_results\$result))

    # R object for other processes to use

    saveRDS(gost_results, file = paste(output_prefix, 'gprofiler2.gost_results.rds', sep = '.'))

    # Write full enrichment table (except parents column as that one throws an error)

    gost_results\$results <- data.frame(
            round_dataframe_columns(gost_results\$result[,-which(names(gost_results\$result) == "parents")], digits=opt\$round_digits),
            check.names = FALSE
        )

    write.table(
        gost_results\$results,
        file = paste(output_prefix, 'gprofiler2', 'all_enriched_pathways_tab', 'tsv', sep = '.'),
        col.names = TRUE,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )

    # Iterate over the enrichment results by source and save separate tables
    for (df in split(gost_results\$result, gost_results\$result\$source)){
        db_source <- df\$source[1]
        df\$short_name <- sapply(df\$term_name, substr, start=1, stop=50)

        df_subset <- data.frame(
            Pathway_name = df\$short_name,
            Pathway_code = df\$term_id,
            DE_genes = df\$intersection_size,
            Pathway_size = df\$term_size,
            Fraction_DE = df\$recall,
            Padj = df\$p_value,
            DE_genes_names = df\$intersection
        )
        df_subset <- data.frame(
            round_dataframe_columns(df_subset, digits=opt\$round_digits),
            check.names = FALSE
        )
        write.table(
            df_subset,
            file = paste(output_prefix, 'gprofiler2', db_source, 'sub_enriched_pathways', 'tsv', sep = '.'),
            col.names = TRUE,
            row.names = FALSE,
            sep = '\t',
            quote = FALSE
        )

        # Enriched pathways horizontal barplots of padj values
        p <- ggplot(df_subset, aes(x=reorder(Pathway_name, Fraction_DE), y=Fraction_DE)) +
            geom_bar(aes(fill=Padj), stat="identity", width = 0.7) +
            geom_text(aes(label=paste0(df_subset\$DE_genes, "/", df_subset\$Pathway_size)), vjust=0.4, hjust=-0.2, size=3) +
            coord_flip() +
            scale_y_continuous(limits = c(0.00, 1.24), breaks = seq(0, 1.24, by = 0.25)) +
            scale_fill_continuous(high = enrich_colors[1], low = enrich_colors[2]) +
            ggtitle(paste("Enriched pathways in", output_prefix)) +
            xlab("") + ylab("Enriched fraction (DE features / Pathway size)")
        ggsave(p, filename = paste(output_prefix, 'gprofiler2', db_source, 'sub_enriched_pathways', 'png', sep = '.'), device = "png", dpi = 300, width=10, limitsize=F)   # Set width to ensure there is enough space for the labels
    }
}

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink("R_sessionInfo.log")
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
gprofiler2.version <- as.character(packageVersion('gprofiler2'))
ggplot2.version <- as.character(packageVersion('ggplot2'))
writeLines(
    c(
        '"\${task.process}":',
        paste('    r-base:', r.version),
        paste('    r-ggplot2:', ggplot2.version),
        paste('    r-gprofiler2:', gprofiler2.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
