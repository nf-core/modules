#!/usr/bin/env Rscript

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
parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Boolean. TRUE if first row is header. False without header.
#' @param row.names The first column is used as row names by default.
#' Otherwise, give another number. Or use NULL when no row.names are present.
#'
#' @return output Data frame
read_delim_flexible <- function(file, header = TRUE, row.names = 1, check.names = TRUE){

    ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1))

    if (ext == "tsv" || ext == "txt") {
        separator <- "\\t"
    } else if (ext == "csv") {
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }

    mat <- read.delim(
        file,
        sep = separator,
        header = header,
        row.names = row.names,
        check.names = check.names
    )
}

#' Get genewise table with logfold changes and connectivity information
#'
#' This function calculates the logfold changes of genes with respect to the reference set,
#' which is dynamically defined as the set of genes that are significantly proportional to
#' each target gene. Note that the output table will only contain genes that are significantly
#' proportional to at least one other gene.
#'
#' @param results Data frame with significant pairs
#' @return Data frame with the following columns:
#'  - lfc = overall logfold change of the gene with respect to the reference set
#'  - lfc_error = median average deviation of the logfold changes -> this reflects the error
#'  - connectivity = size of the reference set -> this also reflects the connectivity of the gene
#   - weighted_connectivity = this reflects the weighted connectivity of the gene, so the lower
#'      the theta the closer to 1 full connectivity. One can also interpret this as the accumulated
#'      between group variance of the gene (as the theta values reflects the between group variance
#'      vs within group variance).
get_genewise_information=function(res,cut=0.01){


    D=ncol(res@counts)
    M=matrix(0,D,4)
    colnames(M)=c("rcDdis","Nsig","lrmD","LFC")
    rownames(M)=colnames(res@counts)
    print("preparing matrices")
    pset=which(res@results[,"FDR"]<cut) #only significant pairs
    re1=prepmat(pset,D, res) #transform pairwise table to DxD matrix #TODO: Ask Ionas if its ok to pass res as an argument

    aset=c(1:nrow(res@results)) #for LFC we need all pairs
    re2=prepmat(aset,D, res)

    print("calculating")
    all=c(1:D)
    for (g in 1:D){
        M[g,"Nsig"]=length(which(res@results[pset,1]==g | res@results[pset,2]==g))

        set1=which(all<g)
        set2=which(all>g)
        M[g,"LFC"]=mean(c(re2[set1,g],-re2[set2,g]))/log(2) #numerator and denominator of lrm swap for genes > g, so lrm2-lrm1

        set=which(re1[all,g]!=0)
        M[g,"lrmD"]=median(c(re1[intersect(set1,set),g],-re1[intersect(set2,set),g]))/log(2) #divide by log2 to get base 2 log

    }

    deg=sort(unique(M[,"Nsig"]),decreasing=TRUE)
    cl=0
    for (i in 1:length(deg)){
        dset=which(M[,"Nsig"]==deg[i])
        cl=cl+length(dset)
        M[dset,"rcDdis"]=cl/D
    }

    ind=order(M[,"rcDdis"])
    M=M[ind,]
    return(M)

}

#this function takes a set of pairs and returns a DxD matrix of lrm1-lrm2:
prepmat=function(subs,D, res){

    pset=subs

    mat1=matrix(0,D,D)
    rownames(mat1)=colnames(res@counts)
    colnames(mat1)=colnames(res@counts)
    mat2=matrix(0,D,D)
    rownames(mat2)=colnames(res@counts)
    colnames(mat2)=colnames(res@counts)
    pa=matrix("",length(pset),2)
    pa[,1]=colnames(res@counts)[res@results[pset,1]]
    pa[,2]=colnames(res@counts)[res@results[pset,2]]

    mat1[pa]=res@results[pset,"lrm1"] #this puts values below the diagonal because of index ordering pa[,1]>pa[,2] in propr
    mat2[pa]=res@results[pset,"lrm2"]
    re=mat1-mat2
    re[upper.tri(re)]=t(re)[upper.tri(re)] #make matrix symmetric
    return(re)
}

#' Plot genewise information
#'
#' This function plots the genewise information, which is a scatter plot of the logfold changes
#' of the genes with respect to the reference set (x-axis) and the accumulated between group variance
#' of the genes (y-axis). The accumulated between group variance is calculated as the sum of 1 - theta
#' values of the genes that are significantly proportional to the target gene. This can be interpreted
#' as the weighted connectivity of the gene in the network.
#'
#' @param results Data frame with genewise information
#' @param output Output png file name
plot_genewise_information <- function(results, output) {

    # Define significant genes
    meanset <- which(results[,"Nsig"] > mean(results[,"Nsig"]))  # Genes above avg significance
    downset <- which(results[,"LFC"] < 0)  # Downregulated genes
    upset <- which(results[,"LFC"] > 0)  # Upregulated genes
    exprob <- max(results[meanset, "rcDdis"])  # Threshold for significance

    # Define tick marks for right-side axis
    Nticks <- rep(0, 4)
    for (t in 0:3) {
        Nticks[t+1] <- max(results[which(-log10(results[,"rcDdis"]) <= t), "Nsig"])
    }

    # Create PNG file
    png(output, height = 12, width = 12, units = "cm", res = 300)

    # Plot empty volcano plot
    plot(
        results[,"LFC"], -log10(results[,"rcDdis"]),
        ylab = "-log10(Probability of Hub)",
        xlab = "Log2 Fold Change",
        type = "n",
        main = "Volcano of Hub Genes"
    )

    # Add reference lines
    abline(v = 0, col = "lightgrey")  # Vertical line at LFC = 0
    abline(h = -log10(exprob), col = "lightgrey")  # Horizontal threshold line

    # Plot all genes in light grey
    points(results[,"LFC"], -log10(results[,"rcDdis"]), pch = 20, cex = 0.5, col = "lightgrey")

    # Highlight significantly downregulated genes (red)
    points(results[intersect(meanset, downset), "LFC"],
           -log10(results[intersect(meanset, downset), "rcDdis"]),
           pch = 20, cex = 0.5, col = "red")

    # Highlight significantly upregulated genes (blue)
    points(results[intersect(meanset, upset), "LFC"],
           -log10(results[intersect(meanset, upset), "rcDdis"]),
           pch = 20, cex = 0.5, col = "blue")

    # Add legend
    legend(
        "bottomright",
        legend = c("Upregulated", "Downregulated", "Below Avg Degree", "Center Gene"),
        col = c("blue", "red", "lightgrey", "black"),
        pch = c(20, 20, 20, 4),
        cex = 0.5
    )

    # Add right-side axis with Nsig values
    axis(side = 4, labels = Nticks, at = c(0:3))

    dev.off()
}


################################################
################################################
## Parse arguments                            ##
################################################
################################################

# Set defaults and classes

opt <- list(
    prefix             = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),

    # input count matrix
    counts             = '$counts',
    features_id_col    = 'gene_id',            # column name of feature ids

    # comparison groups
    samplesheet        = '$samplesheet',
    obs_id_col         = 'sample',             # column name of observation ids
    contrast_variable  = "$contrast_variable", # column name of contrast variable
    reference_group    = "$reference",         # reference group for contrast variable
    target_group       = "$target",            # target group for contrast variable

    # parameters for computing differential proportionality
    alpha              = NA,                   # alpha for boxcox transformation
    moderated          = TRUE,                 # use moderated theta

    # parameters for getting the significant differentially proportional pairs
    fdr                = 0.05,                 # FDR threshold
    permutation        = 0,                    # if permutation > 0, use permutation test to compute FDR
    number_of_cutoffs  = 100,                  # number of cutoffs for permutation test

    # saving options
    # note that pairwise outputs are very large, so it is recommended to save them only when needed
    save_pairwise_full = FALSE,               # save full pairwise results
    save_pairwise      = FALSE,                # save filtered pairwise results
    save_adjacency     = FALSE,                # save adjacency matrix
    save_rdata         = FALSE,                # save rdata

    # other parameters
    seed               = NA,                   # seed for reproducibility
    round_digits       = NA,                   # number of digits to round results
    ncores             = as.integer('$task.cpus')
)

opt_types <- list(
    prefix             = 'character',
    counts             = 'character',
    samplesheet        = 'character',
    features_id_col    = 'character',
    obs_id_col         = 'character',
    contrast_variable  = 'character',
    reference_group    = 'character',
    target_group       = 'character',
    alpha              = 'numeric',
    moderated          = 'logical',
    fdr                = 'numeric',
    permutation        = 'numeric',
    number_of_cutoffs  = 'numeric',
    save_pairwise_full = 'logical',
    save_pairwise      = 'logical',
    save_adjacency     = 'logical',
    save_rdata         = 'logical',
    seed               = 'numeric',
    round_digits       = 'numeric',
    ncores             = 'numeric'
)

# Apply parameter overrides

args_ext <- ifelse('$task.ext.args' == 'null', '', '$task.ext.args')
args_opt <- parse_args(args_ext)
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    } else {

        # Preserve classes from defaults where possible
        args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])

        # handle NA, and avoid errors when NA is provided by user as character
        if (args_opt[[ao]] %in% c('NA', NA)) args_opt[[ao]] <- NA

        # replace values
        opt[[ao]] <- args_opt[[ao]]
    }
}

# Check if required parameters have been provided

required_opts <- c('counts','samplesheet','contrast_variable','reference_group','target_group')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]
if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid

for (file_input in c('counts','samplesheet')){
    if (is.null(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }
    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

# check parameters are valid

if (opt\$permutation < 0) {
    stop('permutation should be a positive integer')
}

print(opt)

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(propr)

################################################
################################################
## Load data                                  ##
################################################
################################################

# set seed when required

if (!is.na(opt\$seed)) {
    warning('Setting seed ', opt\$seed, ' for reproducibility')
    set.seed(opt\$seed)
}

# read input matrix

counts <- read_delim_flexible(
    opt\$counts,
    header = TRUE,
    row.names = opt\$features_id_col,
    check.names = FALSE
)
counts <- t(counts)  # transpose matrix to have features (genes) as columns

# read input samplesheet

samplesheet <- read_delim_flexible(
    opt\$samplesheet,
    header = TRUE,
    row.names = opt\$obs_id_col,
    check.names = FALSE
)

# Check that all samples specified in the input samplesheet are present in the counts
# table. Assuming they are, subset and sort the count table to match the samplesheet

missing_samples <-
    samplesheet[!rownames(samplesheet) %in% rownames(counts), opt\$sample_id_col]

if (length(missing_samples) > 0) {
    stop(paste(
        length(missing_samples),
        'specified samples missing from count table:',
        paste(missing_samples, collapse = ',')
    ))
} else{
    counts <- counts[rownames(samplesheet),] # this will remove non-sample columns, such as metadata columns
    counts <- apply(counts, 2, as.numeric) # if there is a column with non-numeric values, the rest of the matrix will be coerced to character. This will convert it back to numeric
}

# parse group and filter matrix and group values, keeping only the contrasted groups
# TODO propd can also handle more than two groups but that don't work properly with
# the current contrast format. Should we provide an alternative way to do that?

idx <- which(samplesheet[,opt\$contrast_variable] %in% c(opt\$reference_group, opt\$target_group))
counts <- counts[idx,]
samplesheet <- samplesheet[idx,]
group <- as.vector(samplesheet[,opt\$contrast_variable])
group <- as.character(group)
if (length(group) != nrow(counts)) stop('Error when parsing group')
if (length(unique(group)) != 2) stop('Only two groups are allowed for contrast')

################################################
################################################
## Perform differential proportionality       ##
################################################
################################################

# calculate the differential proportionality theta values

pd <- propd(
    counts,
    group    = group,
    alpha    = opt\$alpha,
    weighted = FALSE,
    p        = opt\$permutation
)

# calculate theta moderated, when required
# and calculate F-stat

pd <- updateF(
    pd,
    moderated = opt\$moderated
)
if (opt\$moderated) pd <- setActive(pd, what='theta_mod')

# get significant results based on the FDR-adjusted F-stat p-values, if permutation == 0
# otherwise get them based on the FDR obtained from permutation tests (more computationally expensive but likely more conservative FDRs)

if (opt\$permutation == 0) {

    warning('FDR-adjusted p-values are used to get significant pairs.')

    # get theta value for which FDR is below desired threshold
    # theta_cutoff is FALSE when no theta value has FDR below desired threshold
    # otherwise it is the theta value for which FDR is below desired threshold
    # Only when there is a meaningful theta, we can compute the next steps
    # that involve extracting the significant pairs.

    theta_cutoff <- getCutoffFstat(
        pd,
        pval = opt\$fdr,
        fdr_adjusted = TRUE
    )
    if (theta_cutoff) {

        warning('Significant theta value found: ', theta_cutoff)

        # get adjacency matrix
        # this matrix will have 1s for significant pairs and 0s for the rest
        # diagonals are set to 0

        if (opt\$save_adjacency) {
            adj <- getAdjacencyFstat(
                pd,
                pval = opt\$fdr,
                fdr_adjusted = TRUE
            )
        }

        # get significant pairs

        results_pairwise <- getSignificantResultsFstat(
            pd,
            pval = opt\$fdr,
            fdr_adjusted = TRUE
        )

        # parse genewise information from pairwise results

        results_genewise <- get_genewise_information(pd, cut=0.01) # TODO: Ask Ionas if the value should be provided by the user
    }

} else {

    warning('Permutation tests are used to compute FDR values.')

    # calculate FDR values using permutation tests
    # This test is computationally expensive but it is likely to
    # provide more conservative FDR values.
    # This part will call the updateCutoffs function iteratively
    # as far as it does not find a meaningful theta value
    # and does not reach the maximum number of iterations.

    fdr_table <- data.frame(
        'cutoff' = numeric(0),
        'randcounts' = numeric(0),
        'truecounts' = numeric(0),
        'FDR' = numeric(0)
    )
    theta_cutoff <- FALSE
    max_cutoff <- 1
    ntry <- 0
    while (!theta_cutoff & max_cutoff > 0 & ntry < 10) {
        ntry <- ntry + 1

        # get a list of theta values served as cutoff to calculate the FDR values
        # Given a theta value as cutoff, the FDR is defined as the proportion of
        # false positives obtained from the null distribution vs the total number
        # of positives obtained from the real data.

        cutoffs <- as.numeric(quantile(
            pd@results[pd@results\$theta < max_cutoff, 'theta'],
            seq(0, 1, length.out = opt\$number_of_cutoffs)
        ))

        # update FDR values

        pd <- updateCutoffs(
            pd,
            custom_cutoffs = cutoffs,
            ncores = opt\$ncores
        )
        fdr_table <- rbind(
            pd@fdr[pd@fdr\$cutoff < max_cutoff,],
            fdr_table
        )

        # get theta value for which FDR is below desired threshold
        # theta_cutoff is FALSE when no theta value has FDR below desired threshold
        # otherwise it is the theta value for which FDR is below desired threshold
        # Only when there is a meaningful theta, we can compute the next steps
        # that involve extracting the significant pairs.

        theta_cutoff <- getCutoffFDR(
            pd,
            fdr=opt\$fdr,
            window_size=1
        )

        # update maximum theta value to test the FDR values for the next iteration

        part <- pd@fdr[which(pd@fdr\$truecounts > 0),]
        max_cutoff <- ifelse(nrow(part) > 1, min(part\$cutoff), 0)
    }

    if (theta_cutoff) {

        warning('Significant theta value found: ', theta_cutoff)

        # get adjacency matrix
        # this matrix will have 1s for significant pairs and 0s for the rest
        # diagonals are set to 0

        if (opt\$save_adjacency) {
            adj <- getAdjacencyFDR(
                pd,
                fdr=opt\$fdr,
                window_size=1
            )
        }

        # get significant pairs

        results_pairwise <- getSignificantResultsFDR(
            pd,
            fdr = opt\$fdr,
            window_size = 1
        )

        # parse genewise information from pairwise results

        results_genewise <- get_genewise_information(pd, cut=0.01) # TODO: Ask Ionas if the value should be provided by the user
    }
}

# deal with the situation when no significant thetas are found
# For the moment, we just create empty tables with the same data structure

if (!theta_cutoff) {
    warning('No theta value has FDR below desired threshold.')

    # create empty adjacency matrix

    if (opt\$save_adjacency) {
        adj <- matrix(0, nrow=ncol(counts), ncol=ncol(counts))
        colnames(adj) <- rownames(adj) <- colnames(counts)
    }

    # create empty pairwise results table

    if (opt\$save_pairwise) {
        results <- data.frame(
            'Pair' = character(0),
            'Partner' = character(0),
            'theta' = numeric(0),
            'Fstat' = numeric(0),
            'Pval' = numeric(0),
            'FDR' = numeric(0)
        )
        results_pairwise <- results
    }

    # create empty genewise results table

    results_genewise <- data.frame(
        'features_id_col' = character(0),
        recDdis = numeric(0),
        Nsig = numeric(0),
        lrmD = numeric(0),
        LFC = numeric(0)
    )
    colnames(results_genewise) <- c(opt\$features_id_col, 'recDdis', 'Nsig', 'lrmD', 'LFC')
}

################################################
################################################
## Generate outputs                           ##
################################################
################################################

# save plot of genewise information
# save empty plot if no DE genes were found

if (nrow(results_genewise) > 0) {
    plot_genewise_information(
        results_genewise,
        paste0(opt\$prefix, '.propd.genewise.png')
    )
} else {
    warning('No genewise information to plot.')
    png(paste0(opt\$prefix, '.propd.genewise.png'))
    plot.new()
    dev.off()
}

# save main results - genewise
results_genewise <- as.data.frame(results_genewise)
results_genewise <- results_genewise[order(
    results_genewise\$rcDdis,
    abs(results_genewise\$LFC),
    decreasing = TRUE
),]

if (!is.na(opt\$round_digits)) {
    cols <- sapply(results_genewise, is.numeric)
    results_genewise[,cols] <- round(
        results_genewise[,cols],
        digits = opt\$round_digits
    )
}

write.table(
    results_genewise,
    file      = paste0(opt\$prefix, '.propd.genewise.tsv'),
    col.names = TRUE,
    row.names = FALSE,
    sep       = '\\t',
    quote     = FALSE
)

# save rdata, if required

if (opt\$save_rdata) {
    saveRDS(
        pd,
        file = paste0(opt\$prefix, '.propd.rds')
    )
}

# save pairwise results, if required

if (opt\$save_pairwise) {

    # unfiltered pairwise results table

    if (opt\$save_pairwise_full) {
        results <- getResults(pd)
        rm(pd)
        results <- results[order(
            results\$theta,
            results\$FDR
        ), c('Pair', 'Partner', 'theta', 'Fstat', 'Pval', 'FDR')]

        if (!is.na(opt\$round_digits)) {
            cols <- sapply(results, is.numeric)
            results[,cols] <- round(
                results[,cols],
                digits = opt\$round_digits
            )
        }

        write.table(
            results,
            file      = paste0(opt\$prefix, '.propd.pairwise.tsv'),
            col.names = TRUE,
            row.names = FALSE,
            sep       = '\\t',
            quote     = FALSE
        )
    }

    # filtered pairwise results table

    results_pairwise <- results_pairwise[order(
        results_pairwise\$theta,
        results_pairwise\$FDR
    ), c('Pair', 'Partner', 'theta', 'Fstat', 'Pval', 'FDR')]

    if (!is.na(opt\$round_digits)) {
        cols <- sapply(results_pairwise, is.numeric)
        results_pairwise[,cols] <- round(
            results_pairwise[,cols],
            digits = opt\$round_digits
        )
    }

    write.table(
        results_pairwise,
        file      = paste0(opt\$prefix, '.propd.pairwise_filtered.tsv'),
        col.names = TRUE,
        row.names = FALSE,
        sep       = '\\t',
        quote     = FALSE
    )
}

# save adjacency matrix, if required

if (opt\$save_adjacency) {
        write.table(
            adj,
            file      = paste0(opt\$prefix, '.propd.adjacency.csv'),
            col.names = TRUE,
            row.names = TRUE,
            sep       = ',',
            quote     = FALSE
        )
}

# save FDR values, if permutation tests were run

if (opt\$permutation > 0) {
    fdr_table <- fdr_table[order(fdr_table\$cutoff),]

    if (!is.na(opt\$round_digits)) {
        fdr_table\$FDR <- round(
            fdr_table\$FDR,
            digits = opt\$round_digits
        )
    }

    write.table(
        fdr_table,
        file      = paste0(opt\$prefix, '.propd.fdr.tsv'),
        col.names = TRUE,
        row.names = FALSE,
        sep       = '\\t',
        quote     = FALSE
    )
}

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste0(opt\$prefix, ".R_sessionInfo.log"))
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

propr.version <- as.character(packageVersion('propr'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-propr:', propr.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
