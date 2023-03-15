process COUNTRNA_COLLECT {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each groupData

    output:
        tuple val(groups), path("FRASER_output"),   emit: result

    shell:
        groups      = groupData.group
        FR_output   = groupData.output
        grSplit     = groupData.grSplit
        grNonSplit  = groupData.grNonSplit
        spliceSites = groupData.spliceSites

        rawCountsJ  = groupData.rawCountsJ
        rawCountsSS = groupData.rawCountsSS
        '''
            #!/usr/bin/env Rscript --vanilla

            suppressPackageStartupMessages({
                library(rmarkdown)
                library(knitr)
                library(devtools)
                library(yaml)
                library(BBmisc)
                library(GenomicAlignments)
                library(tidyr)
                library(data.table)
                library(dplyr)
                library(plotly)
                library(DelayedMatrixStats)
                library(FRASER)
                library(rhdf5)
                library(purrr)
                library(base)
            })

            dataset    <- "!{groups}"
            saveDir    <- "FRASER_output"

            # Read FRASER object
            file.copy("!{FR_output}", ".", recursive = TRUE)
            fds <- loadFraserDataSet(dir="FRASER_output", name=paste0("raw-local-", dataset))
            splitCounts_gRanges <- readRDS("!{grSplit}")
            spliceSiteCoords <- readRDS("!{spliceSites}")

            # Get splitReads and nonSplitRead counts in order to store them in FRASER object
            splitCounts_h5 <- HDF5Array::HDF5Array("!{rawCountsJ}", "rawCountsJ")
            splitCounts_se <- SummarizedExperiment(
            colData = colData(fds),
            rowRanges = splitCounts_gRanges,
            assays = list(rawCountsJ=splitCounts_h5)
            )


            nonSplitCounts_h5 <- HDF5Array::HDF5Array("!{rawCountsSS}", "rawCountsSS")
            nonSplitCounts_se <- SummarizedExperiment(
            colData = colData(fds),
            rowRanges = spliceSiteCoords,
            assays = list(rawCountsSS=nonSplitCounts_h5)
            )

            # Add Counts to FRASER object
            fds <- addCountsToFraserDataSet(fds=fds, splitCounts=splitCounts_se,
                                            nonSplitCounts=nonSplitCounts_se)

            # Save final FRASER object
            fds <- saveFraserDataSet(fds)
        '''
}
