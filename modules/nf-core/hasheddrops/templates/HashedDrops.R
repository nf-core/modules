#!/usr/bin/env Rscript

################################################
################################################
## Fucntions                                  ##
################################################
################################################

# Helper function for NULL condition
string_to_null <- function(x, val = "NULL") if (x == val) NULL else x
null_to_string <- function(x, val = "NULL") if (is.null(x)) val else x

string_to_logical <- function(input) {
  if (input == "FALSE") {
    FALSE
  } else if (input == "TRUE") {
    TRUE
  } else {
    stop(paste0(input, " is not a valid logical. Use 'FALSE' or 'TRUE'."))
  }
}

################################################
################################################
## USE PARAMETERS FROM NEXTFLOW               ##
################################################
################################################

# cast parameters from nextflow

# hashedDrops parameters
hto_matrix <- '$hto_matrix'
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
library(DropletUtils) # for hashedDrops() and emptyDrops()

################################################
################################################
## Main Process                               ##
################################################
################################################

hto <- Read10X(data.dir = hto_matrix, gene.column = gene_col)

# die bekomme ich mit runempty Drops
#is.cell <- NULL

# determine hto_input and ambient_input
if (runEmptyDrops) {

    rna <- Read10X(data.dir = rna_matrix, gene.column = gene_col)

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

    # which droplets are actual cells
    is.cell <- emptyDrops_out\$FDR <= isCellFDR
    hto_input <- hto[, which(is.cell)]

    if (ambient) {
        ambient_input <- metadata(emptyDrops_out)\$ambient
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
  confident.min = confidentMin,
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
  "confidentMin",
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
  confidentMin,
  null_to_string(combinations)
)

params <- data.frame(Argument, Value)
write.csv(params, paste0(prefix ,"_params_hasheddrops.csv"))

#--------- save emptyDrops() results  ---------#

# create a plot with results from emptyDrops() or save an empty png
png(paste0(prefix, "_emptyDrops.png"))
if(runEmptyDrops){
    colors <- ifelse(is.cell, "red", "black")
    plot(emptyDrops_out\$Total, -emptyDrops_out\$LogProb, col = colors, xlab = "Total UMI count", ylab = "-Log Probability")
}else{
    plot.new()

}
dev.off()

write.csv(emptyDrops_out,paste0(prefix, "_emptyDrops.csv"))
saveRDS(emptyDrops_out,file = paste0(prefix, "_emptyDrops.rds"))

#--------- save hashedDrops() results ---------#

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
    col = colors,
    xlab = "Log-fold change between best and second HTO",
    ylab = "Log-fold change between second HTO and ambient"
    )
}else{

    plot.new()
}

dev.off()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
seurat.version <- as.character(packageVersion('Seurat'))
dropletutils.version <- as.character(packageVersion('DropletUtils'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    r-seurat:', seurat.version),
        paste('    dropletutils:', dropletutils.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
