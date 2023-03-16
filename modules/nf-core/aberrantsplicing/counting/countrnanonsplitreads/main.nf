process COUNTRNA_NONSPLITREADS_SAMPLEWISE {
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
        tuple val(grData), path("FRASER_output")    , emit: result

    shell:
        groups       = groupData.group
        sampleID     = groupData.sampleID
        grSplit      = groupData.grSplit
        FR_output    = groupData.output

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
        # Read FRASER object
        file.copy("!{FR_output}", ".", recursive = TRUE)
        fds <- loadFraserDataSet(dir="FRASER_output", name=paste0("raw-local-", dataset))

        # Get sample id from wildcard
        sample_id <- "!{sampleID}"


        # Read splice site coordinates from RDS
        spliceSiteCoords <- readRDS("!{grSplit}")

        # Count nonSplitReads for given sample id
        sample_result <- countNonSplicedReads(sample_id,
                                            splitCountRanges = NULL,
                                            fds = fds,
                                            minAnchor=5,
                                            recount="!{recount}",
                                            spliceSiteCoords=spliceSiteCoords,
                                            longRead="!{longRead}")

        message(date(), ": ", dataset, ", ", sample_id,
                " no. splice junctions (non split counts) = ", length(sample_result))

        # run the version part
        cat(file="versions.yml", "!{task.process}:\naberrantexpression: 1.3.0")
        '''
}
