process COUNTRNA_NONSPLITREADS_MERGE {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each groupData
        val recount
        val longRead

    output:
        tuple val(grData), path("FRASER_output"), \
            path("FRASER_output/savedObjects/raw-local-*/rawCountsSS.h5")    , emit: result

    shell:
        groups       = groupData.group
        FR_output    = groupData.output.join(",")
        grNonSplit   = groupData.grNonSplit

        grData = groupData

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

        options("FRASER.maxSamplesNoHDF5"=1)
        options("FRASER.maxJunctionsNoHDF5"=-1)

        h5disableFileLocking()

        dataset    <- "!{groups}"
        recount    <- "!{recount}" == "true"

        # copy all folders together
        for (x in strsplit("!{FR_output}", ",")[[1]]) {
            file.copy(x, ".", recursive=TRUE, overwrite=TRUE)
        }

        # Read FRASER object
        fds <- loadFraserDataSet(dir="FRASER_output", name=paste0("raw-local-", dataset))

        # Read splice site coordinates from RDS
        splitCounts_gRanges <- readRDS("!{grNonSplit}")

        # If samples are recounted, remove the merged ones
        nonSplitCountsDir <- file.path("FRASER_output", "savedObjects",
                                paste0("raw-local-", dataset), 'nonSplitCounts')
        if(recount == TRUE & dir.exists(nonSplitCountsDir)){
        unlink(nonSplitCountsDir, recursive = TRUE)
        }

        # Get and merge nonSplitReads for all sample ids
        nonSplitCounts <- getNonSplitReadCountsForAllSamples(fds=fds,
                                                            NcpuPerSample=2,
                                                            splitCountRanges=splitCounts_gRanges,
                                                            minAnchor=5,
                                                            recount=FALSE,
                                                            longRead="!{longRead}")

        message(date(), ":", dataset, " nonSplit counts done")

        # run the version part
        cat(file="versions.yml", "!{task.process}:\naberrantexpression: 1.3.0")
        '''
}
