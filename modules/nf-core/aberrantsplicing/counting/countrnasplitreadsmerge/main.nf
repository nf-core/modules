process COUNTRNA_SPLITREADS_MERGE {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each sampleIDs
        val recount
        val minExpressionInOneSample

    output:
        tuple val(grData), path("FRASER_output"), path("FRASER_output/cache/raw-local-*/gRanges_splitCounts.rds"), \
            path("FRASER_output/cache/raw-local-*/gRanges_NonSplitCounts.rds"), \
            path("FRASER_output/cache/raw-local-*/spliceSites_splitCounts.rds"),
            path("FRASER_output/savedObjects/raw-local-*/rawCountsJ.h5")  , emit: result

    shell:
        groups = sampleIDs.group
        output = sampleIDs.output.join(",")

        grData = sampleIDs

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

            # copy all folders together
            for (x in strsplit("!{output}", ",")[[1]]) {
                file.copy(x, ".", recursive=TRUE, overwrite=TRUE)
            }

            dataset    <- "!{groups}"
            recount    <- "!{recount}" == "true"
            minExpressionInOneSample <- !{minExpressionInOneSample}

            # Read FRASER object
            fds <- loadFraserDataSet(dir="FRASER_output", name=paste0("raw-local-", dataset))

            # If samples are recounted, remove the merged ones
            splitCountsDir <- file.path("FRASER_output", "savedObjects",
                                        paste0("raw-local-", dataset), 'splitCounts')
            if(recount == TRUE & dir.exists(splitCountsDir)){
                unlink(splitCountsDir, recursive = TRUE)
            }

            # Get and merge splitReads for all sample ids
            splitCounts <- getSplitReadCountsForAllSamples(fds=fds,
                                                        NcpuPerSample=2,
                                                        recount=FALSE)
            # Extract, annotate and save granges
            splitCountRanges <- rowRanges(splitCounts)

            # Annotate granges from the split counts
            dir.create("FRASER_output/cache/raw-local-!{groups}")

            splitCountRanges <- FRASER:::annotateSpliceSite(splitCountRanges)
            saveRDS(splitCountRanges, "FRASER_output/cache/raw-local-!{groups}/gRanges_splitCounts.rds")

            # Create ranges for non split counts
            # Subset by minExpression
            maxCount <- rowMaxs(assay(splitCounts, "rawCountsJ"))
            passed <- maxCount >= minExpressionInOneSample
            # extract granges after filtering
            splitCountRanges <- splitCountRanges[passed,]

            saveRDS(splitCountRanges, "FRASER_output/cache/raw-local-!{groups}/gRanges_NonSplitCounts.rds")

            # Extract splitSiteCoodinates: extract donor and acceptor sites
            # take either filtered or full fds
            spliceSiteCoords <- FRASER:::extractSpliceSiteCoordinates(splitCountRanges, fds)
            saveRDS(spliceSiteCoords, "FRASER_output/cache/raw-local-!{groups}/spliceSites_splitCounts.rds")

            message(date(), ": ", dataset, " total no. splice junctions = ",
                    length(splitCounts))

            # run the version part
            cat(file="versions.yml", "!{task.process}:\naberrantexpression: 1.3.0")
        '''
}
