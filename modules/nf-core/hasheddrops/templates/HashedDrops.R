#!/usr/bin/env Rscript

################################################
################################################
## Fucntions                                  ##
################################################
################################################

# Helper function for NULL condition
string_to_null <- function(x, val = "NULL") if (x == val) NULL else x
null_to_string <- function(x, val = "NULL") if (is.null(x)) val else x

<<<<<<< HEAD
string_to_logical <- function(input) {
  if (input == "FALSE") {
    FALSE
  } else if (input == "TRUE") {
    TRUE
  } else {
    stop(paste0(input, " is not a valid logical. Use 'FALSE' or 'TRUE'."))
  }
}

=======
>>>>>>> 11d1bc9b0 (save changes)
################################################
################################################
## USE PARAMETERS FROM NEXTFLOW               ##
################################################
################################################

# cast parameters from nextflow

# hashedDrops parameters
hto_matrix <- '$hto_matrix'
<<<<<<< HEAD
rna_matrix <- '$rna_matrix'
lower <- as.numeric('$lower')
niters <- as.numeric('$niters')
testAmbient <- string_to_logical('$testAmbient')
ignore <- string_to_null('$ignore')
alpha <- string_to_null('$alpha')
round <- string_to_logical('$round')
byRank <- string_to_null('$byRank')
isCellFDR <- as.numeric('$isCellFDR')
ambient <- string_to_logical('$ambient')
minProp <- as.numeric('$minProp')
pseudoCount <- as.numeric('$pseudoCount')
constantAmbient <- string_to_logical('$constantAmbient')
doubletNmads <- as.numeric('$doubletNmads')
doubletMin <- as.numeric('$doubletMin')
doubletMixture <- string_to_logical('$doubletMixture')
confidentNmads <- as.numeric('$confidentNmads')
confidentMin <- as.numeric('$confidentMin')
combinations <- string_to_null('$combinations')
runEmptyDrops <- string_to_logical('$runEmptyDrops')
gene_col <- as.numeric('$gene_col')
prefix <- '$prefix'

print("------------------------")
print('$runEmptyDrops')
print(runEmptyDrops)
print("------------------------")

=======
hto_matrix <- '$rna_matrix'
lower <- as.numeric('$lower')
niters <- as.numeric('$niters')
testAmbient <- as.logical('$testAmbient')
ignore_hashedDrops <- string_to_null('$ignore_hashedDrops')
alpha_hashedDrops <- string_to_null('$alpha_hashedDrops')
round <- as.logical('$round')
byRank <- string_to_null('$byRank')
isCellFDR <- as.numeric('$isCellFDR')
ambient <- as.logical('$ambient')
minProp <- as.numeric('$minProp')
pseudoCount <- as.numeric('$pseudoCount')
constantAmbient <- as.logical('$constantAmbient')
doubletNmads <- as.numeric('$doubletNmads')
doubletMin <- as.numeric('$doubletMin')
doubletMixture <- as.logical('$doubletMixture')
confidentNmads <- as.numeric('$confidentNmads')
confidentMin <- as.numeric('$confidentMin')
combinations <- string_to_null('$combinations')
runEmptyDrops <-  as.logical('$runEmptyDrops')
gene_col <- as.numeric('$gene_col')
prefix <- '$prefix'

>>>>>>> 11d1bc9b0 (save changes)
# check if the file exists
if (! file.exists(hto_matrix)){
    stop(paste0(hto_matrix, ' is not a valid file'))
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(Seurat)  # for Read10X()
<<<<<<< HEAD
library(DropletUtils) # for hashedDrops() and emptyDrops()
=======
library(cellhashR) # for hashedDrops() and emptyDrops()
>>>>>>> 11d1bc9b0 (save changes)

################################################
################################################
## Main Process                               ##
################################################
################################################

hto <- Read10X(data.dir = hto_matrix, gene.column = gene_col)

# die bekomme ich mit runempty Drops
#is.cell <- NULL

# determine hto_input and ambient_input
<<<<<<< HEAD
print("some output:: ")
print(runEmptyDrops)
=======
>>>>>>> 11d1bc9b0 (save changes)
if (runEmptyDrops) {

    rna <- Read10X(data.dir = rna_matrix, gene.column = gene_col)

<<<<<<< HEAD
    #print(rna)

=======
>>>>>>> 11d1bc9b0 (save changes)
    emptyDrops_out <- emptyDrops(
    rna,
    lower = lower,
    niters = niters,
    test.ambient = testAmbient,
    ignore = NULL,
    alpha = alpha,
    round = round,
    by.rank = byRank
    )

<<<<<<< HEAD
    #print(emptyDrops_out)

    # which droplets are actual cells
    is.cell <- emptyDrops_out\$FDR <= isCellFDR
    hto_input <- hto[, which(is.cell)]

    if (ambient) {
        ambient_input <- metadata(emptyDrops_out)\$ambient
=======
    # which droplets are actual cells
    is.cell <- emptyDrops_out$FDR <= isCellFDR
    hto_input <- hto[, which(is.cell)]

    if (ambient == TRUE) {
        ambient_input <- metadata(emptyDrops_out)$ambient
>>>>>>> 11d1bc9b0 (save changes)
    } else {
        ambient_input <- NULL
    }
} else {
    ambient_input <- NULL
    hto_input <- hto

    # only important for saving the results
    emptyDrops_out <- data.frame()
}

hashedDrops_out <- hashedDrops(
  hto_input,
  min.prop = minProp,
  ambient = ambient_input,
  pseudo.count = pseudoCount,
  constant.ambient = constantAmbient,
  doublet.nmads = doubletNmads,
  doublet.min = doubletMin,
  doublet.mixture = doubletMixture,
  confident.nmads = confidentNmads,
<<<<<<< HEAD
  confident.min = confidentMin,
=======
  confident.min = confidenMin,
>>>>>>> 11d1bc9b0 (save changes)
  combinations = combinations
)

################################################
################################################
## SAVING RESULTS                             ##
################################################
################################################

#----- saving parameters in a dataframe  ------#
Argument <- c(
  "hto_matrix",
  "lower",
  "niters",
  "testAmbient",
  "ignore",
  "alpha",
  "round",
  "byRank",
  "isCellFDR",
  "ambient",
  "minProp",
  "pseudoCount",
  "constantAmbient",
  "doubletNmads",
  "doubletMin",
  "doubletMixture",
  "confidentNmads",
<<<<<<< HEAD
  "confidentMin",
=======
  "confidenMin",
>>>>>>> 11d1bc9b0 (save changes)
  "combinations"
)

Value <- c(
  hto_matrix,
  lower,
  niters,
  testAmbient,
  null_to_string(ignore),
  null_to_string(alpha),
  round,
  null_to_string(byRank),
  isCellFDR,
  ambient,
  null_to_string(minProp),
  pseudoCount,
  constantAmbient,
  doubletNmads,
  doubletMin,
  doubletMixture,
  confidentNmads,
<<<<<<< HEAD
  confidentMin,
=======
  confidenMin,
>>>>>>> 11d1bc9b0 (save changes)
  null_to_string(combinations)
)

params <- data.frame(Argument, Value)
write.csv(params, paste0(prefix ,"_params_hasheddrops.csv"))

#--------- save emptyDrops() results  ---------#

# create a plot with results from emptyDrops() or save an empty png
png(paste0(prefix, "_emptyDrops.png"))
if(runEmptyDrops){
    colors <- ifelse(is.cell, "red", "black")
<<<<<<< HEAD
    plot(emptyDrops_out\$Total, -emptyDrops_out\$LogProb, col = colors, xlab = "Total UMI count", ylab = "-Log Probability")
=======
    plot(emptyDrops_out$Total, -emptyDrops_out$LogProb, col = colors, xlab = "Total UMI count", ylab = "-Log Probability")
>>>>>>> 11d1bc9b0 (save changes)
}else{
    plot.new()

}
dev.off()

write.csv(emptyDrops_out,paste0(prefix, "_emptyDrops.csv"))
saveRDS(emptyDrops_out,file = paste0(prefix, "_emptyDrops.rds"))

#--------- save hashedDrops() results ---------#

<<<<<<< HEAD
write.csv(params, paste0(prefix, "_params_hasheddrops.csv"))
write.csv(hashedDrops_out,paste0(prefix,"_results_hasheddrops.csv"))
saveRDS(hashedDrops_out,file = paste0(prefix,"_hasheddrops.rds"))

png(paste0(prefix, "_plot_hasheddrops.png"))
if (sum(is.na(hashedDrops_out\$LogFC2)) != length(hashedDrops_out\$LogFC2)) {

    colors <- ifelse(hashedDrops_out\$Confident,
    "black",
    ifelse(hashedDrops_out\$Doublet, "red", "grey")
    )

    plot(
    hashedDrops_out\$LogFC,
    hashedDrops_out\$LogFC2,
=======
write.csv(params, paste0(prefix, "params_hasheddrops.csv"))
write.csv(hashedDrops_out,paste0(prefix,"results_hasheddrops.csv"))
saveRDS(hashedDrops_out,file = paste0(prefix,"hasheddrops.rds"))

png(paste0(args$outputdir, "/", "plot_hasheddrops.png"))
if (sum(is.na(hashedDrops_out$LogFC2)) != length(hashedDrops_out$LogFC2)) {

    colors <- ifelse(hashedDrops_out$Confident,
    "black",
    ifelse(hashedDrops_out$Doublet, "red", "grey")
    )

    plot(
    hashedDrops_out$LogFC,
    hashedDrops_out$LogFC2,
>>>>>>> 11d1bc9b0 (save changes)
    col = colors,
    xlab = "Log-fold change between best and second HTO",
    ylab = "Log-fold change between second HTO and ambient"
    )
<<<<<<< HEAD
}else{

    plot.new()
}

dev.off()

=======

    dev.off()
}

>>>>>>> 11d1bc9b0 (save changes)
################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
seurat.version <- as.character(packageVersion('Seurat'))
<<<<<<< HEAD
dropletutils.version <- as.character(packageVersion('DropletUtils'))
=======
cellhashR.version <- as.character(packageVersion('cellhashR'))
>>>>>>> 11d1bc9b0 (save changes)

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    r-seurat:', seurat.version),
<<<<<<< HEAD
        paste('    dropletutils:', dropletutils.version)
=======
        paste('    cellhashR:', cellhashR.version)
>>>>>>> 11d1bc9b0 (save changes)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
