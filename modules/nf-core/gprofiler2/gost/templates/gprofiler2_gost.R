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
        columns <- colnames(df)[(unlist(lapply(df, is.numeric), use.names=F))]    # extract only numeric columns for rounding
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
    de_file = '$de_file',
    de_id_column = 'gene_id',
    organism = NULL,
    sources = NULL,
    output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    significant = T,
    measure_underrepresentation = F,
    correction_method = 'gSCS',
    evcodes = F,
    pval_threshold = 0.05,
    gmt_file = '$gmt_file',
    token = NULL,
    background_file = '$background_file',
    background_column = NULL,
    domain_scope = 'annotated',
    min_diff = 1,
    round_digits = -1,
    palette_name = 'Blues',
    archive = 'gprofiler'
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
required_opts <- c('output_prefix')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]

if (length(missing) > 0) {
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}
if (is.null(opt\$organism) && opt\$gmt_file == "" && is.null(opt\$token)) {
    stop('Please provide organism, gmt_file or token.')
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

output_prefix <- paste0(opt\$output_prefix, ".gprofiler2")

# Create empty output table as it is a mandatory output
file.create(paste(output_prefix, 'all_enriched_pathways', 'tsv', sep = '.'))

if (nrow(de.genes) > 0) {

    query <- de.genes[[opt\$de_id_column]]

    ################################################
    ################################################
    # Run gprofiler processes and generate outputs #
    ################################################
    ################################################

    set.seed(1) # This will ensure that reruns have the same plot colors

    sources <- opt\$sources
    if (!is.null(sources)) {
        sources <-  strsplit(opt\$sources, split = ",")[[1]]
    }
    if (!is.null(sources)) {
        sources <-  strsplit(opt\$sources, split = ",")[[1]]
    }

    if (!is.null(opt\$token)) {

        # First check if a token was provided
        token <- opt\$token

    } else if (!is.null(opt\$organism)) {

        # Next, check if organism was provided. Get the GMT file from gprofiler and save both the full file as well as the filtered one to metadata
        base_url <- paste0("https://biit.cs.ut.ee/", opt\$archive)
        gmt_url <- paste0(base_url, "//static/gprofiler_full_", opt\$organism, ".ENSG.gmt")
        set_base_url(base_url)
        tryCatch(
            {
                gmt_path <- paste0("gprofiler_full_", opt\$organism, ".ENSG.gmt")
                if (!is.null(sources)) {
                    gmt_path <- paste0("gprofiler_full_", opt\$organism, ".", paste(sources, collapse="_"), ".ENSG_filtered.gmt")
                }
                download <- download.file(gmt_url, gmt_path)
                if (download != 0) {
                    print("Failed to fetch the GMT file from gprofiler with this URL:")
                    print(gmt_url)
                    print("For reproducibility reasons, try to download the GMT file manually by visiting https://biit.cs.ut.ee/gprofiler/gost, then selecting the correct organism and, in datasources, clicking 'combined ENSG.gmt'.")
                } else {
                    if (!is.null(sources)) {
                        gmt <- Filter(function(line) any(startsWith(line, sources)), readLines(gmt_path))
                        print(paste0("GMT file successfully downloaded and filtered. Please note that for some sources, the GMT file may not contain any entries as these cannot be retrieved from gprofiler; in this case, the GMT file may be completely empty."))
                        writeLines(gmt, gmt_path)
                    }
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
        token <- opt\$organism

    } else {

        # Last option: Use custom GMT file
        gmt_path <- opt\$gmt_file

        # If sources are set, extract only requested entries (gprofiler will NOT filter automatically!)
        if (!is.null(sources)) {
            gmt <- Filter(function(line) any(startsWith(line, sources)), readLines(opt\$gmt))
            gmt_path <- paste0(strsplit(basename(opt\$gmt_file), split = "\\\\.")[[1]][[1]], ".", paste(sources, collapse="_"), "_filtered.gmt")
            writeLines(gmt, gmt_path)
        }
        token <- upload_GMT_file(gmt_path)

        # Add gost ID to output GMT name so that it can be reused in future runs
        file.rename(gmt_path, paste0(strsplit(basename(opt\$gmt_file), split = "\\\\.")[[1]][[1]], ".", paste(sources, collapse="_"), "_gostID_", token, "_filtered.gmt"))

    }


    # If custom background_file was provided, read it
    if (opt\$background_file != "") {
        intensities_table <- read_delim_flexible(
            file = opt\$background_file
        )
        # If only 1 col, it is a list, not a matrix
        if (ncol(intensities_table) == 1) {
            background <- intensities_table[,1]                                 # Extract first column from df
            background <- append(background, colnames(intensities_table)[1])    # First entry was put into header, add it to vector
        } else {
            # Otherwise it's a matrix
            # Set rownames to background_column if param was set
            if (!is.null(opt\$background_column)) {
                if (opt\$background_column %in% colnames(intensities_table)) {
                    rownames(intensities_table) <- intensities_table[[opt\$background_column]]
                    intensities_table[[opt\$background_column]] <- NULL
                } else {
                    stop(paste0("Invalid background_column argument: ", opt\$background_column,
                                ". Valid columns are: ", paste(colnames(intensities_table), collapse=", "), "."))
                }
            } else {

            # Otherwise set rownames to first column
                rownames(intensities_table) <- intensities_table[,1]
                intensities_table <- intensities_table[,-1]
            }

            # Rownames are set, now remove non-numeric columns
            nums <- unlist(lapply(intensities_table, is.numeric), use.names = FALSE)
            intensities_table <- intensities_table[, nums]
            # Keep only rownames which have abundance
            background <- rownames(subset(intensities_table, rowSums(intensities_table, na.rm = TRUE)>0))
        }
    } else {
        background <- NULL
    }

    # Name the query as it will otherwise be called 'query_1' which will also determine the gostplot title
    q <- list(query)
    names(q) <- c(output_prefix)
    gost_results <- gost(
        query=q,
        organism=token,
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
                    file = paste(output_prefix, 'gostplot', 'html', sep = '.')
                    )

        # Create a static plot and save to PNG
        static_plot <- gostplot(gost_results, capped=T, interactive=F)
        ggsave(plot = static_plot, filename = paste(output_prefix, 'gostplot', 'png', sep = '.'), width = 10, height = 7)

        # Subset gost results to those pathways with a min. number of differential features
        gost_results\$result <- gost_results\$result[which(gost_results\$result\$intersection_size>=opt\$min_diff),]

        # annotate query size (number of differential features in contrast)
        gost_results\$result\$original_query_size <- rep(length(as.character(de.genes\$Ensembl_ID)), nrow(gost_results\$result))

        # R object for other processes to use

        saveRDS(gost_results, file = paste(output_prefix, 'gost_results.rds', sep = '.'))

        # Write full enrichment table (except parents column as that one throws an error)

        gost_results\$results <- data.frame(
                round_dataframe_columns(gost_results\$result[,-which(names(gost_results\$result) == "parents")], digits=opt\$round_digits),
                check.names = FALSE
            )

        write.table(
            gost_results\$results,
            file = paste(output_prefix, 'all_enriched_pathways', 'tsv', sep = '.'),
            col.names = TRUE,
            row.names = FALSE,
            sep = '\t',
            quote = FALSE
        )

        # Iterate over the enrichment results by source and save separate tables
        for (df in split(gost_results\$result, gost_results\$result\$source)){

            db_source <- df\$source[1]
            df_subset <- data.frame(
                Pathway_name = df\$term_name,
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
                file = paste(output_prefix, db_source, 'sub_enriched_pathways', 'tsv', sep = '.'),
                col.names = TRUE,
                row.names = FALSE,
                sep = '\t',
                quote = FALSE
            )

            # For plot, shorten pathway names as they can get quite long (full name can be looked up in the table)
            df_subset\$Pathway_name <- sapply(df_subset\$Pathway_name, substr, start=1, stop=50)

            # Extract 3 colors from the chosen palette (2 are sufficient, but brewer.pal has a minimum of 3); first and last will be used for plot
            colors <- RColorBrewer::brewer.pal(3, opt\$palette_name)

            # Enriched pathways horizontal barplots of padj values
            p <- ggplot(df_subset, aes(x=reorder(Pathway_name, Fraction_DE), y=Fraction_DE)) +
                geom_bar(aes(fill=Padj), stat="identity", width = 0.7) +
                geom_text(aes(label=paste0(df_subset\$DE_genes, "/", df_subset\$Pathway_size)), vjust=0.4, hjust=-0.2, size=3) +
                theme(plot.title.position = "plot") +
                coord_flip() +
                scale_y_continuous(limits = c(0.00, 1.24), breaks = seq(0, 1.24, by = 0.25)) +
                ggtitle(paste("Enriched", db_source, "pathways")) +
                xlab("") + ylab("Enriched fraction (DE features / Pathway size)") +
                scale_fill_continuous(high = colors[1], low = colors[3])

            # Save plot with set width to ensure there is enough space for the labels; adapt height to nrow but limit it to 100 as there will be an error for too high values
            ggsave(p, filename = paste(output_prefix, db_source, 'sub_enriched_pathways', 'png', sep = '.'), device = "png", width=10, height=min(100, 1.5+nrow(df_subset)*0.15), limitsize=F)
        }
    }
} else {
    print("No differential features found, pathway enrichment analysis with gprofiler2 will be skipped.")
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
        '"$task.process":',
        paste('    r-base:', r.version),
        paste('    r-ggplot2:', ggplot2.version),
        paste('    r-gprofiler2:', gprofiler2.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
