process FILTER {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each groupData
        val sampleAnnotation
        val params

    output:
        tuple val(groups), path("FRASER_output"), emit: result

    shell:
        groups       = groupData.group
        FR_output    = groupData.output

        minDeltaPsi  = params.minDeltaPsi
        minExpressionInOneSample = params.minExpressionInOneSample
        filter = params.filter

        exCountIDs   = groupData.exIDs.join(",")
        exCountFiles = groupData.exFiles.join(",")
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

            opts_chunk$set(fig.width=12, fig.height=8)

            # copy FRASER_output
            file.copy("!{FR_output}", ".", recursive = TRUE)

            # input
            workingDir <- "FRASER_output"
            dataset    <- "!{groups}"
            exCountIDs <- strsplit("!{exCountIDs}", ",")[[1]]
            exCountFiles <- strsplit("!{exCountFiles}", ",")[[1]]
            sample_anno_file <- "!{sampleAnnotation}"
            minExpressionInOneSample <- !{minExpressionInOneSample}
            minDeltaPsi <- !{minDeltaPsi}

            exCountIDs <- strsplit("!{exCountIDs}", ",")[[1]]
            exCountFiles <- strsplit("!{exCountFiles}", ",")[[1]]

            fds <- loadFraserDataSet(dir="FRASER_output", name=paste0("raw-local-", dataset))

            # Add external data if provided by dataset
            if(length(exCountIDs) > 0){
                message("create new merged fraser object")
                fds <- saveFraserDataSet(fds,dir = workingDir, name=paste0("raw-", dataset))

                for(resource in unique(exCountFiles)){
                    exSampleIDs <- exCountIDs[exCountFiles == resource]
                    exAnno <- fread(sample_anno_file, key="RNA_ID")[J(exSampleIDs)]
                    setnames(exAnno, "RNA_ID", "sampleID")

                    ctsNames <- c("k_j", "k_theta", "n_psi3", "n_psi5", "n_theta")
                    ctsFiles <- paste0(resource, "/", ctsNames, "_counts.tsv.gz")

                    # Merging external counts restricts the junctions to those that
                    # are only present in both the counted (fromBam) junctions AND the
                    # junctions from the external counts.
                    fds <- mergeExternalData(fds=fds, countFiles=ctsFiles,
                            sampleIDs=exSampleIDs, annotation=exAnno)
                    fds@colData$isExternal <- as.factor(!is.na(fds@colData$SPLICE_COUNTS_DIR))
                }
            } else {
                message("symLink fraser dir")
                file.symlink(paste0(workingDir, "savedObjects/","raw-local-", dataset),
                            paste0(workingDir, "savedObjects/","raw-", dataset))

                fds@colData$isExternal <- as.factor(FALSE)
                workingDir(fds) <- workingDir
                name(fds) <- paste0("raw-", dataset)
            }

            # filter for expression and write it out to disc.
            fds <- filterExpressionAndVariability(fds,
                    minExpressionInOneSample = minExpressionInOneSample,
                    minDeltaPsi = minDeltaPsi,
                    filter=FALSE)

            devNull <- saveFraserDataSet(fds,dir = workingDir)

            # Keep junctions that pass filter
            name(fds) <- dataset
            if ("!{filter}" == "true") {
                filtered <- mcols(fds, type="j")[,"passed"]
                fds <- fds[filtered,]
                message(paste("filtered to", nrow(fds), "junctions"))
            }

            seqlevels(fds) <- seqlevelsInUse(fds)
            colData(fds)$sampleID <- as.character(colData(fds)$sampleID)
            fds <- saveFraserDataSet(fds,dir = workingDir)

            # run the version part
            cat(file="versions.yml", "!{task.process}:\naberrantexpression: 1.3.0")
        '''
}
