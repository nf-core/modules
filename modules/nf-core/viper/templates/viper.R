#!/usr/bin/env Rscript

library(viper)

# Inject variables
exp_mat_file <- "${expression_matrix}"
network_file <- "${network}"

# Load regulon and expression matrix
exp_mat <- as.matrix(read.table(exp_mat_file, header=TRUE))
regul <- aracne2regulon(network_file, exp_mat)

# Run VIPER
vpres <- viper(exp_mat, regul, verbose = FALSE)

# Output result
output_file <- file.path("viper_results.tsv")
write.table(vpres, file = output_file, sep = "\t", quote = FALSE)

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

viper_version <- as.character(packageVersion('viper'))

writeLines(
    c(
        '"${task.process}":',
        paste('    viper:', viper_version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################