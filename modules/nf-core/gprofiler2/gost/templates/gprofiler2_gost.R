#!/usr/bin/env Rscript

# Written by Oskar Wacker (https://github.com/WackerO) in
# collaboration with Stefan Czemmel (https://github.com/qbicStefanC)
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
    de_table = '$de_table',
    de_id_column = 'gene_id'
    contrast_variable = '$contrast_variable',
    reference_level = '$reference',
    target_level = '$target',
    blocking_variables = NULL,
    genome = '$genome',
    significant = T,
    min_diff = 1,
    correction_method = 'fdr',
    sources = NULL,
    evcodes = F,
    user_threshold = '0.05',
    background = NULL,
    gmt = NULL,
    counts_table = NULL,
    domain_scope = 'annotated',
    round_digits = -1
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

required_opts <- c('contrast_variable', 'reference_level', 'target_level', 'output_prefix', 'sources')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]

if (length(missing) > 0) {
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid

for (file_input in c('de_table')) {
    if (is.null(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }

    if (! file.exists(opt[[file_input]])) {
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

# Determine organism and libary

org_names <- list(
    'GRCh37' = 'hsapiens',
    'GRCh38' = 'hsapiens',
    'GRCm38' = 'mmusculus',
    'TAIR10' = 'athaliana'
)
org_keytypes <- list(
    'GRCh37' = 'ENSEMBL',
    'GRCh38' = 'ENSEMBL',
    'GRCm38' = 'ENSEMBL',
    'TAIR10' = 'TAIR'
)
org_libraries <- list(
    'GRCh37' = 'org.Hs.eg.db',
    'GRCh38' = 'org.Hs.eg.db',
    'GRCm38' = 'org.Mm.eg.db',
    'TAIR10' = 'org.At.tair.db'
)

org_name <- org_names[[opt\$genome]]
org_keytype <- org_keytypes[[opt\$genome]]
org_library <- org_libraries[[opt\$genome]]

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

#library(plotly)
library(gprofiler2)

#BiocManager::install(org_library, version="3.17", force=T)
#library(org_library, character.only=T)
#species_library_installed <- get(opt\$species_library)

################################################
################################################
## READ IN DIFFERENTIAL GENES FILE            ##
################################################
################################################

de.genes <-
    read_delim_flexible(
        file = opt\$de_table
    )
if (nrow(de.genes) == 0) {
  opt <- options(show.error.messages = FALSE) # Do not throw error so that rest of pipeline continues
  print("No differential features found, pathway enrichment analysis with gprofiler2 will be skipped.")
  stop()
}
query <- de.genes[[de_id_column]]
# TODO IS AS.CHARACTER NEEDED ABOVE?


################################################
################################################
# Run gprofiler processes and generate outputs #
################################################
################################################

output_prefix <- opt\$output_prefix

# Generate plots for all requested normalizations; also, save normalized protein groups for limma















if (isProvided(opt\$custom_gmt)){

    # If custom GMT file was provided, extract only requested datasources (gprofiler will NOT filter automatically!)
    gmt <- (Filter(function(line) any(startsWith(line, datasources)), readLines(opt\$gmt)))
    gmt_file <- paste0(strsplit(basename(opt\$gmt), split = "\\\\.")[[1]], 1), "_filtered.gmt"))
    writeLines(out_gmt, gmt_file)
    gost_id <- upload_GMT_file(out_path)
} else {

    # Otherwise, get the GMT file from gprofiler and save both the full file as well as the filtered one to metadata
    gmt_url <- paste0("https://biit.cs.ut.ee/gprofiler//static/gprofiler_full_", organism, ".ENSG.gmt")
    tryCatch(
        {
            wget_command <- paste0("wget ", gmt_url)
            sys_return <- system(wget_command)
            if (sys_return != 0 && !("saved" %in% sys_return)) {
                print("Failed to fetch the GMT file from gprofiler with this URL:")
                print(gmt_url)
                print("For reproducibility reasons, try to download the GMT file manually by visiting https://biit.cs.ut.ee/gprofiler/gost, then selecting the correct organism and, in datasources, clicking 'combined ENSG.gmt'.")
            }
        },
        error=function(gost_error) {
            print("Failed to fetch the GMT file from gprofiler with this URL:")
            print(gmt_url)
            print("Got error:")
            print(gost_error)
            print("For reproducibility reasons, please try to download the GMT file manually by visiting https://biit.cs.ut.ee/gprofiler/gost, then selecting the correct organism and, in datasources, clicking 'combined ENSG.gmt'. Then provide it to the pipeline with the parameter `--custom_gmt`")
        }
    )
    gost_id <- organism
}




















# If custom background was provided, use that
if (!is.null(opt\$background)) {
    background <- readLines(opt\$background)
}
# Else if counts table was provided, use as background all features whose row sums are > 0
else if (!is.null(opt\$counts_table)) {
    counts_table <- read_delim_flexible(
        file = opt\$counts_table,
        row.names = 1
    )
    background <- rownames(counts_table)[rowSums(counts(cds))>0]
}

gostres <- gost(
    query=query,
    organism=gost_id,
    significant=opt\$significant,
    correction_method=opt\$correction_method,
    sources=opt\$sources,
    evcodes=opt\$evcodes,
    user_threshold=opt\$user_threshold,
    custom_bg=background,
    domain_scope=opt\$domain_scope
)

if (nrow(gostres$result) > 0) {

    # Create interactive plot and save to HTML
    interactive_plot <- gostplot(gostres, capped=T, interactive=T)
    interactive_plot[['x']][['layout']][['annotations']][[1]][['x']] <- -opt\$adj_pval_threshold

    # limit gostplot y maximum dynamically for all subplots
    #for (counter in c(1:length(contrast_files))) {

    #    expression <- parse(text=paste0("pg2 %>% layout(yaxis", ifelse(counter>1, counter, "")," = list(range = list(0, y_max)))"))
    #    pg2 <- eval(expression)
    #}

    #add export button to plot
    #plot <- config(plot, modeBarButtonsToAdd = list(svg_exp)) %>% layout(width=500, height = 400*length(contrast_files))

    # Save interactive plot as HTML
    htmlwidgets::saveWidget(
                widget = interactive_plot,
                file = file.path(output_prefix, "gostplot.html")
                )

    # Create a static plot and save to PNG

    static_plot <- gostplot(gostres, capped=T, interactive=F)
    png(filename = file.path(output_prefix, "gostplot.png"), width = 900, height = 600)
    static_plot
    dev.off()

    # ggsave(
    #     static_plot,
    #     filename = "gostplot.png",
    #     device="png",
    #     height=10,
    #     width=15,
    #     units="cm",
    #     dpi=300,
    #     limitsize=F
    #     )

} else { 
    opt <- options(show.error.messages = FALSE) # Do not throw error so that rest of pipeline continues
    print("No enriched features found by gprofiler2.")
    stop()
}




# Subset gost results to those pathways with a min. number of differential features
gostres$result <- gostres$result[which(gostres$result$intersection_size>=opt/$min_diff),]

# annotate query size (number of differential features in contrast)
gostres$result$original_query_size <- rep(length(as.character(DE_genes$Ensembl_ID)), nrow(gostres$result))



# Iterate over the enrichment results by source
for (df in split(pathway_gostres, pathway_gostres$source)){
    db_source <- df$source[1]
    df$short_name <- sapply(df$term_name, substr, start=1, stop=50)

    df_subset <- data.frame(
        Pathway_name = df$short_name,
        Pathway_code = df$term_id,
        DE_genes = df$intersection_size,
        Pathway_size = df$term_size,
        Fraction_DE = (df$intersection_size / df$term_size),
        Padj = df$p_value,
        DE_genes_names = df$intersection
    )
    df_subset <- data.frame(
        round_dataframe_columns(df_subset, digits=opt\$round_digits),
        check.names = FALSE
    )
    write.table(
        out_df,
        file = paste(output_prefix, 'gprofiler2', db_source, 'enriched_pathways_tab', 'tsv', sep = '.'),
        col.names = TRUE,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )

    # Enriched pathways horizontal barplots of padj values
    p <- ggplot(df_subset, aes(x=reorder(Pathway_name, Fraction_DE), y=Fraction_DE)) +
        geom_bar(aes(fill=Padj), stat="identity", width = 0.7) +
        geom_text(aes(label=paste0(df_subset$DE_genes, "/", df_subset$Pathway_size)), vjust=0.4, hjust=-0.5, size=3) +
        coord_flip() +
        scale_y_continuous(limits = c(0.00, 1.00)) +
        scale_fill_continuous(high = "#132B43", low = "#56B1F7") +
        ggtitle("Enriched pathways") +
        xlab("") + ylab("Enriched fraction (DE features / Pathway size)")
    ggsave(p, filename = paste(output_prefix, 'gprofiler2', db_source, 'enriched_pathways_tab', 'png', sep = '.'), device = "png", dpi = 300, limitsize=F)   #height = 5+0.5*nrow(df_subset), units = "cm", 
}




























org_name



# R object for other processes to use

saveRDS(proteinGroups, file = paste(output_prefix, 'proteus.raw_proteingroups.rds', sep = '.'))

# Remove parents column to be able to save the table in tsv format

#gostres$result$parents <- NULL

# Write enrichment table without parents column (otherwise will throw error)

out_df <- data.frame(
        round_dataframe_columns(gostres$result[,names(gostres$result) != "parents"], digits=opt\$round_digits),
        check.names = FALSE
    )
out_df[[opt\$protein_id_col]] <- rownames(proteinGroups\$tab) # proteus saves the IDs as rownames; save these to a separate column
out_df <- out_df[c(opt\$protein_id_col, colnames(out_df)[colnames(out_df) != opt\$protein_id_col])] # move ID column to first position


write.table(
    out_df,
    file = paste(output_prefix, 'gprofiler2', 'enriched_pathways_tab', 'tsv', sep = '.'),
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

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
limma.version <- as.character(packageVersion('limma'))
plotly.version <- as.character(packageVersion('plotly'))
proteus.version <- as.character(packageVersion('proteus'))
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    r-proteus-bartongroup:', proteus.version),
        paste('    r-plotly:', plotly.version),
        paste('    bioconductor-limma:', limma.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
