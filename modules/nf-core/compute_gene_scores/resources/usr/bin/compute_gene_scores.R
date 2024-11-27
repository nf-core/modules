#!/usr/bin/env Rscript

#Libraries
library(singscore)
library(GSVA)
library(stringr)
library(dplyr)
library(yaml)
library(FactoMineR)
library(stringr)
library(tibble)

###Gene Set Functions

#  helper functions
loadMappingTable <- function(ref.file="./ensembl_86_mpens.txt") {
    s2e <- read.delim(file=ref.file,
                    stringsAsFactors = F, header=F, skip=1)
    colnames(s2e) <- c("ensembl", "genesymbol")
    return(s2e)
}

geneSymbol2Ensembl <- function(x, s2e=NULL) {
    if(is.null(s2e)) s2e <- loadMappingTable()
    # also check for correct colnames (TODO)
    ens <- unique(s2e[s2e$genesymbol %in% x,"ensembl"])
    # print(ens)
    if(length(ens)<length(x)) warning(paste0(length(x)-length(ens), " symbols could not be mapped."))
    return(ens)
}

geneSet2Ensembl <- function(sets, s2e=NULL) {
    sets.ensembl <- lapply(sets, geneSymbol2Ensembl, s2e)
    return(sets.ensembl)
}

geneSet2List <- function(set) {
    set.list <- list()
    set.list$upSet <- set
    set.list$downSet <- NULL
    return(set.list)
}

geneSets2Lists <- function(sets) {
    sets.tmp <- lapply(sets, geneSet2List)
    return(sets.tmp)
}

flattenGenesets <- function(set) {
    if(is.vector(set)) return(set)
    if(is.null(set$downSet)) return(set[['upSet']])
    else(return(unique(c(set[['upSet']], set[['downSet']]))))
}

# flatten gene set for undirected algorithms
geneSetAsVector <- function(set) {
    if(is.character(set)) {
        return(set)
    } else if(is.data.frame(set)) {
        return(set[,1])
    } else {
        stop("Gene Set must be a character vector or a data.frame")
    }
}

# Create a pseudo-network as input for PAI algorithm
geneSetAsNetwork <- function(set) {
    if(is.character(set)) {
        set <- set
    } else if(is.data.frame(set)) {
        set <- set[,1]
    } else {
        stop("Gene Set must be a character vector or a data.frame")
    }
    set.net <- data.frame(id1=set, id2=set, stringsAsFactors = F)
    return(set.net)
}

# Create a list for input to Singscore calculations
geneSetAsList <- function(set) {
    retSet <- list(upSet=character(),
                    downSet=character(),
                    knownDir=FALSE)
    if(is.character(set)) {
        retSet$upSet <- set
        retSet$downSet <- NULL
        kownDir <- FALSE
    } else if(is.data.frame(set)) {
    if(!is.numeric(set[,2])) {
        stop("Second column of gene set data frame must be numeric.")
        } else if(all(set[,2]==0)) {
        retSet$upSet <- set[,1]
        retSet$downSet <- NULL
        kownDir <- FALSE
    }
        else {
        retSet$upSet <- set[set[,2]>0,1]
        retSet$downSet <- set[set[,2]<0,1]
        retSet$knownDir <- TRUE
    }
    } else {
        stop("Gene Set must be a character vector or a data.frame")
    }
    if(length(retSet$upSet)==0) retSet$upSet<-NULL
    if(length(retSet$downSet)==0) retSet$downSet<-NULL
    return(retSet)
}

#Scores
## Create an actual Geneset object from GSEABase package
## not used at the moment
createGeneset <- function(set, name.set) {
    gs <- GeneSet(unique(set), setName=name.set)
    return(gs)
}


##
Singscore <- function(d, x, ...) {
    if(is.list(x) && is.null(x$downSet)) s <- simpleScore(d, upSet = x$upSet, knownDirection=x$knownDir, ...)
    if(is.list(x) && !is.null(x$downSet)) s <- simpleScore(d, upSet = x$upSet, downSet=x$downSet, knownDirection = x$knownDir, ...)
    return(s)
}


## Calculate singscores (use log TPMs or RPKMs)
calculateSingscores <- function(xp, genes, raw=FALSE, ...) {

    geneSets <- lapply(genes, geneSetAsList)

    ranked.exp <- rankGenes(xp)

    sscore <- lapply(geneSets, function(x) Singscore(ranked.exp, x, ...))

    flattenSingScores <- function(x,n) {
    if("TotalScore" %in% colnames(x)) {
        one.row <- x[,"TotalScore",drop=F]
        return (one.row)
    }
    }

    sscore.df <- lapply(names(sscore), function(x) flattenSingScores(sscore[[x]],x))
    names(sscore.df) <- names(sscore)
    sscore.f <- bind_cols(sscore.df)
    colnames(sscore.f) <- names(sscore.df)
    rownames(sscore.f) <- rownames(sscore.df[[1]])
    sscore.f <- t(sscore.f)
    if(raw) return(sscore)
    else return(sscore.f)
}

## Calculate GSVA (check kcdf paramater dependent on data)
calculateGSVAscores <- function(xp, genes, ...) {

    genes.flat <- lapply(genes, geneSetAsVector)

    if(is.data.frame(xp)) xp.matrix <- as.matrix.data.frame(xp)
    if(is.matrix(xp)) xp.matrix <- xp

    scores <- gsva(xp.matrix, genes.flat, method="gsva", ...)
    return(scores)
}


## Calculate SSGSEA (use log TPMs or RPKMs)
calculateSSGSEAscores <- function(xp, genes, ...) {

    genes.flat <- lapply(genes, geneSetAsVector)

    if(is.data.frame(xp)) xp.matrix <- as.matrix.data.frame(xp)
    if(is.matrix(xp)) xp.matrix <- xp

    scores <- gsva(xp.matrix, genes.flat, method="ssgsea", ...)
    return(scores)
}

## Calculate PLAGE (probably also TPMs or RPKMs)
calculatePLAGEscores <- function(xp, genes, correct.dir=FALSE, ...) {

    genes.flat <- lapply(genes, geneSetAsVector)

    if(is.data.frame(xp)) xp.matrix <- as.matrix.data.frame(xp)
    if(is.matrix(xp)) xp.matrix <- xp

    plagescores <- gsva(xp.matrix, genes.flat, method="plage", ...)

    if(correct.dir) {
        meanscores <- calculateMeanScores(xp.matrix, genes.flat)
        i <- intersect(rownames(meanscores), rownames(plagescores))
    if(length(i)!=length(rownames(plagescores))) stop("Correction failed")
        plagescores <- plagescores*sign(sapply(i, function(x) cor(plagescores[x,],meanscores[x,])))
    }
    return(plagescores)
}


calculateMedianScores <- function(xp, genes, ...) {

    genes.flat <- lapply(genes, geneSetAsVector)

    if(is.data.frame(xp)) xp.matrix <- as.matrix.data.frame(xp)
    if(is.matrix(xp)) xp.matrix <- xp

    median.result <- t(sapply(genes.flat, function(x) apply(xp.matrix[x, ],2,median, ...)))
    return(median.result)
}

calculateMeanScores <- function(xp, genes, ...) {

    genes.flat <- lapply(genes, geneSetAsVector)

    if(is.data.frame(xp)) xp.matrix <- as.matrix.data.frame(xp)
    if(is.matrix(xp)) xp.matrix <- xp

    mean.result <- t(sapply(genes.flat, function(x) apply(xp.matrix[x, ],2,mean, ...)))
    return(mean.result)
}




## geometrix mean functions from stackoverflow
## https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
## Note, normalized read counts can be < 1 => log2 cpm can be  < 0 => geometric mean of log2 cpm values can be infinite and is not well defined
## We implement that anyway, because it is described in Herazo-Maya et al 2017, a publication from the Kaminski lab:
## https://pubmed.ncbi.nlm.nih.gov/28942086/

# Geometric mean
# Naive implementation assuming all values are > 0.
gm_mean_naive = function(a){prod(a)^(1/length(a))}

# Geometric mean
# Implementation including only values > 0.
gm_mean_simple = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Geometric mean
# Implementation, which allows options.
gm_mean_advanced = function(x, na.rm=TRUE, zero.propagate = FALSE){
    if(any(x < 0, na.rm = TRUE)){
        return(NaN)
    }
    if(zero.propagate){
        if(any(x == 0, na.rm = TRUE)){
        return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
    } else {
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
}

# wrapper called by calculateSamsscore to adapt up and down regulation information
Samsscore <- function(xp, x, center_rows = TRUE,geom_mean   = FALSE, ...) {
    if(is.list(x) && !is.null(x$downSet))  s <- sams_score(xp, up_genes = x$upSet,
                                                        down_genes = x$downSet,
                                                        center_rows = center_rows,
                                                        geom_mean   = geom_mean,
                                                        ...)
    if(is.list(x) && is.null(x$downSet))   s <- sams_score(xp, up_genes = x$upSet,
                                                        down_genes = character(length = 0L),
                                                        center_rows = center_rows,
                                                        geom_mean   = geom_mean,
                                                        ...)
return(s)
}

# sams score wrapper to match implementation of other scores
calculateSamsscores <- function(xp, genes, raw = FALSE,
                                which_col = "score_sum",
                                center_rows = TRUE,
                                geom_mean   = FALSE,
                                ...) {

geneSets <- lapply(genes, geneSetAsList)

samsscore <- lapply(geneSets, function(x) Samsscore(xp, x,
                                                    center_rows = center_rows,
                                                    geom_mean   = geom_mean,
                                                    ...))
if (raw) return(samsscore)

my_flat <- function(x) {
    my_tib <- dplyr::tibble(sample = rownames(samsscore[[x]]),
                            !!x   := samsscore[[x]][, which_col])
    my_tib
}

score_flat <- lapply(names(samsscore), my_flat)

score_flat_join <- Reduce(function(...) dplyr::full_join(..., by = "sample"), score_flat)

ret_val <- score_flat_join %>% data.table::as.data.table() %>%
    data.table::melt(id.vars = "sample", variable.name = "geneset", value.name = "score") %>%
    data.table::dcast(geneset~sample, value.var = "score") %>%
    as.matrix(rownames = "geneset")
ret_val[, colnames(xp), drop = FALSE]
}


# basic implementation of samsscore
sams_score <- function(expr_mat,
                    up_genes = character(),
                    down_genes = character(),
                    center_rows = TRUE,
                    geom_mean   = FALSE,
                    method = "sams") { # method "sams" or "median", median is just the column median

    expr_mat <- as.matrix(expr_mat)
    assertthat::assert_that(is.numeric(expr_mat), msg = "expression matix is not numeric")

    assertthat::assert_that(all(up_genes %in% rownames(expr_mat)),
                            msg = "At least one element in up_genes is not in rownames of expression matrix.")
    assertthat::assert_that(all(down_genes %in% rownames(expr_mat)),
                            msg = "At least one element down_genes is not in rownames of expression matrix.")
    if (center_rows) {
        if (geom_mean) {
            row_means <- apply(expr_mat, 1, gm_mean_simple)
        } else {
            row_means <- rowMeans(expr_mat, na.rm = TRUE)
        }
        expr_mat  <- sweep(expr_mat, 1, row_means, "-")
    }

    ## select data for signature
    up_matrix   <- expr_mat[up_genes, ]
    down_matrix <- expr_mat[down_genes, ]

    ##set values, which do not follow up and down-reglation expection to 0
    up_cens   <- up_matrix
    up_cens[up_matrix < 0]     <- matrix(0, nrow = nrow(up_matrix), ncol = ncol(up_matrix))[up_matrix < 0]

    down_cens <- down_matrix
    down_cens[down_matrix > 0] <- matrix(0, nrow = nrow(down_matrix), ncol = ncol(down_matrix))[down_matrix > 0]

    ##calculate fraction of up and down-regulated genes
    #up_fraction   <- apply(up_matrix,   2, function(x) sum(x >= 0) / length(x))
    up_fraction   <- apply(up_matrix,   2, function(x) sum(x > 0) / length(x))
    down_fraction <- apply(down_matrix, 2, function(x) sum(x < 0) / length(x))

    up_score   <- colSums(up_cens)   * up_fraction
    down_score <- colSums(down_cens) * down_fraction

    if(method == "median") {
        up_score   <- apply(up_matrix, 2, median, na.rm = TRUE)
        down_score <- apply(down_matrix, 2, median, na.rm = TRUE)
    }

    my_score_mat <- cbind(up_score = up_score, down_score = down_score)

    my_row_sums  <- rowSums(abs(my_score_mat), na.rm = TRUE)

    if(method == "median") {
        my_score_mat_mod <- cbind(up_score = up_score, down_score = down_score * (-1))
        my_row_sums  <- rowSums(my_score_mat_mod, na.rm = TRUE)
    }

    both_missing <- apply(my_score_mat, 1, function(x) all(is.na(x)))
    my_row_sums[both_missing] <- NA

    my_score_mat <- cbind(my_score_mat, score_sum = my_row_sums)

    return(my_score_mat)
}


## check signature

checkSignatures <- function(xp, genes) {
    g.v <- lapply(genes, geneSetAsVector)
    names(g.v) <- names(genes)

    genes.xp <- unique(rownames(xp))
    genes.diff <- lapply(g.v, setdiff, genes.xp)
    genes.inter <- lapply(g.v, intersect, genes.xp)
    genes.length <- unlist(lapply(genes.diff, length))

    genes.diff <- lapply(genes.diff, paste0, collapse=", ")
    genes.diff <- unlist(genes.diff)

    genes.inter <- lapply(genes.inter, paste0, collapse=", ")
    genes.inter <- unlist(genes.inter)

    l <- unlist(lapply(g.v, length))

    ret.df <- tibble(geneset=names(g.v),
                    size.geneset=l,
                    size.missing=genes.length,
                    size.left=l-genes.length,
                    genes.missing=genes.diff,
                    genes.used=genes.inter)
    return(ret.df)
}


## signature overlap warning
lengthWarning <- function(x) {
    warning("Geneset: ", x[1], ": ", x[3], " of ", x[2], " genes missing: ", x[5], call.=F)
}

## global wrapper function

calculateSignatureScores <- function(xp, genes, method, phenotypes=NULL, case=NULL, control=NULL, raw=FALSE, ...) {
    # Report mismatches between matrix and signatures
    checks <- checkSignatures(xp, genes)
    # print(checks)
    # TODO: new function to send warning

    checks.filtered <- checks %>% filter(size.missing>0)
    apply(checks.filtered, 1, lengthWarning)


    gene.sets.out <- checks %>% filter(size.left==0) %>% pull(geneset)
    if(length(gene.sets.out)>0) warning("The following genesets show no overlap with the expression matrix: " , paste0(gene.sets.out,collapse=", "), call.=F)

    # Keep only genesets with at least 1 gene overlapping
    gene.sets.left <- checks %>% filter(size.left>0) %>% pull(geneset)
    genes <- genes[gene.sets.left]

    # Stop if now genesets are left
    if(length(genes)==0) stop("No genesets show an overlap with the expression matrix. Please check if identifiers match.")

    # Use intersection of gene set with expression matrix for further calculations
    genes <- lapply(genes, intersect, rownames(xp))


# Run actual score calculations
switch(method,
        plage = calculatePLAGEscores(xp, genes, ...),
        plage.dir = calculatePLAGEscores(xp, genes, correct.dir = TRUE, ...),
        gsva = calculateGSVAscores(xp, genes, ...),
        singscore = calculateSingscores(xp, genes, raw, ...),
        ssgsea = calculateSSGSEAscores(xp, genes, ...),
        median = calculateMedianScores(xp, genes, ...),
        mean = calculateMeanScores(xp, genes, ...),
        sams = calculateSamsscores(xp, genes, ...),
        stop(method, " is not a valid method for score calculation."))
}


## wrapper for multiple methods

multipleSignatureScores <- function(xp, genes, methods) {
    stopifnot((is.data.frame(xp) || is.matrix(xp)),
            is.list(genes),
            length(genes)>0,
            is.character(methods), length(methods)>0)

    scores <- lapply(methods, function(x) calculateSignatureScores(xp, genes, x))
    names(scores) <- methods
    return(scores)
}



## Check consistency across scores

# add suffix to rownames
reformatRownames <- function(x, n){
    rownames(x) <- rownames(x) %>% str_replace_all("\\s|\\-",".") %>% paste0(".",n)
    return(x)
}


# generic function to generate a correlation matrix, preferably based on output of multipleSignatureScores
correlationMatrix <- function(score.list, names.list=NULL) {
    stopifnot(is.list(score.list))

    if(is.null(names.list)) names.list <- names(score.list)
    stopifnot(length(score.list)==length(names.list))
    renamed.score.list <- list()
    for(i in 1:length(names.list)) {
        renamed.score.list[[i]] <- reformatRownames(score.list[[i]], names.list[i])
    }
    names(renamed.score.list) <- names(score.list)

    combined.matrix <- do.call(rbind, renamed.score.list)

    cor.matrix.all <- cor(t(combined.matrix))

    return(cor.matrix.all)
}

# legacy function
correlationMatrixOld <- function(plageScores = NULL, gsvaScores = NULL, ssgseaScores = NULL, paiScores = NULL, singScores = NULL) {

    if(!is.null(plageScores)) plageScores <- reformatRownames(plageScores,"PLAGE")
    if(!is.null(gsvaScores)) gsvaScores  <- reformatRownames(gsvaScores, "GSVA")
    if(!is.null(ssgseaScores)) ssgseaScores <- reformatRownames(ssgseaScores, "SSGSEA")
    if(!is.null(paiScores)) paiScores <- reformatRownames(paiScores, "PAI")
    if(!is.null(singScores))  singScores <- reformatRownames(singScores, "Singscore")

    combined.matrix <- rbind(plageScores, gsvaScores, ssgseaScores, paiScores, singScores)

    cor.matrix.all <- cor(t(combined.matrix))

    return(cor.matrix.all)
}


## Comparisons between groups

## Boxplots

boxplotScores <- function(score.matrix, set, group, title="") {
    bp.df <- data.frame(score=score.matrix[set,], group=group)
    ggplot(bp.df, aes(x=group, y=score)) +
        geom_boxplot() +
        geom_point() +
        ggtitle(title)
}


sessionInfo()


####Commandline Argument parsing###
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("Usage: process_nanostring_data.R <scores.yaml> <counts.tsv> algorithm", call.=FALSE)
}
####Uncomment if debugging

gene_sets_to_compute <- args[1]
counts_input <- args[2]
score.method <- args[3]

#Input parsing yaml
yaml_parsed <- read_yaml(gene_sets_to_compute)

#Get GEX into shape for Nanostring pipeline
counts <- read.table(counts_input, sep="\t", check.names=F, header=T) %>%
    filter(`CodeClass` == "Endogenous") %>%
    select(-`CodeClass`) %>%
    column_to_rownames("Name")

# explicit signature check against expression matrix to generate an QC overview
cs <- checkSignatures(counts, yaml_parsed)

#Write to multiQC output

qc.file="signature_scores_qc_mqc.txt"
line="# id: nf-core-nanoflow-signature-score-qc
# section_name: 'Signature Score QC'
# description: 'Compare Signatures to Expression Matrix'
# plot_type: 'table'
# section_href: 'https://github.com/nf-core/nanoflow'"
write(line,
    file = qc.file)
write.table(cs,
            file = qc.file,
            append=TRUE,
            sep="\t",
            row.names = FALSE,
            quote = FALSE,
            na=""
            )

## Compute scores we need in our case
scores <- calculateSignatureScores(xp = counts, genes = yaml_parsed, method=score.method)
s.r <- 1
s.c <- 1
if(ncol(scores)>1) {
    hc <- hclust(dist(t(scores)))
    s.c <- hc$order
}
if(nrow(scores)>1) {
    hc.r <- hclust(dist(scores))
    s.r <- hc.r$order
}
scores.df <- as.data.frame.matrix(t(scores[s.r,s.c]))
scores.df <- rownames_to_column(scores.df, "sample")

#Write to multiQC output
score.file="signature_scores_mqc.txt"
line=paste0("# id: nf-core-nanoflow-signature-score
# section_name: 'Signature Scores'
# description: 'Signature Scores: Algorithm: ",score.method,"'
# plot_type: 'heatmap'
# section_href: 'https://github.com/nf-core/nanoflow'")
write(line, file = score.file)

write.table(scores.df,
            file = score.file,
            append=TRUE,
            sep="\t",
            row.names = FALSE,
            quote = FALSE,
            na=""
)



