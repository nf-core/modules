process MERGECOUNTS {
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each (groupData)
        path (sampleAnn)

    output:
        tuple val(preprocess), val(groupName), path ("total_counts.Rds")    , emit: result
        path "versions.yml"                                                 , emit: versions

    shell:
        groupName = groupData.group

        counts = groupData.rds.join(",")
        count_ranges = groupData.preprocess.countRanges

        rnaIDs = groupData.rnaID.join(",")

        preprocess = groupData.preprocess
        '''
            #!/usr/bin/env Rscript --vanilla

            suppressPackageStartupMessages({
                library(data.table)
                library(dplyr)
                library(BiocParallel)
                library(SummarizedExperiment)
            })
            count_files <- strsplit("!{counts}", ",")[[1]]
            register(MulticoreParam(1))
            count_ranges <- readRDS("!{count_ranges}")

            # Read counts
            selectedIDs <- strsplit("!{rnaIDs}", ",")[[1]]
            counts_list <- bplapply(count_files, function(f){
                if(grepl('Rds$', f))
                    assay(readRDS(f))
                else {
                    ex_counts <- as.matrix(fread(f), rownames = "geneID")

                    index = which(count_files == f)[1]

                    # remove the existing ID from a given index
                    selectedID <- selectedIDs[index]
                    selectedIDs <<- selectedIDs[-index]

                    subset(ex_counts, select = selectedID)
                }
            })

            # check rownames and proceed only if they are the same
            row_names_objects <- lapply(counts_list, rownames)
            if( length(unique(row_names_objects)) > 1 ){
            stop('The rows (genes) of the count matrices to be merged are not the same.')
            }

            # merge counts
            merged_assays <- do.call(cbind, counts_list)
            total_counts <- SummarizedExperiment(assays=list(counts=merged_assays))
            colnames(total_counts) <- gsub('.bam', '', colnames(total_counts))

            # assign ranges
            rowRanges(total_counts) <- count_ranges

            # Add sample annotation data (colData)
            sample_anno <- fread("!{sampleAnn}",
                                colClasses = c(RNA_ID = 'character', DNA_ID = 'character'))
            sample_anno <- sample_anno[, .SD[1], by = RNA_ID]
            col_data <- data.table(RNA_ID = as.character(colnames(total_counts)))
            col_data <- left_join(col_data, sample_anno, by = "RNA_ID")
            rownames(col_data) <- col_data$RNA_ID
            colData(total_counts) <- as(col_data, "DataFrame")
            rownames(colData(total_counts)) <- colData(total_counts)$RNA_ID

            # save in RDS format
            saveRDS(total_counts, "total_counts.Rds")

            # run the version part
            cat(file="versions.yml", "!{task.process}:\naberrantexpression: 1.3.0")
        '''
}
