process DEFINE_DATASET_FROM_ANNO {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each group
        val bams
        path (sampleAnnotation)

    output:
        tuple val(group), path("*.tsv")           , emit: result

    shell:
        bamsStr = bams.join(",")
        group = group
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

            ## helper functions
            write_tsv <- function(x, file, row.names = FALSE, ...){
                write.table(x=x, file=file, quote=FALSE, sep='\t', row.names= row.names, ...)
            }

            splits <- strsplit("!{bamsStr}", ",")[[1]]
            IDs <- splits[seq(1, length(splits), 2)]
            bamFile <- splits[seq(2, length(splits), 2)]

            mapping <- data.frame(IDs, bamFile)
            names(mapping) <- c("sampleID", "bamFile")

            #+ input
            name          <- "!{group}"
            annoFile      <- "!{sampleAnnotation}"

            #+ dataset name

            anno    <- fread(annoFile)

            subset_ids <- IDs
            annoSub <- anno[RNA_ID %in% subset_ids]
            annoSub <- subset(annoSub, grepl(name, DROP_GROUP))
            setnames(annoSub, "RNA_ID", "sampleID")
            colData <- merge(annoSub,
                mapping)
            setcolorder(colData, unique(c("sampleID", "STRAND", "PAIRED_END", "bamFile"), colnames(annoSub)))

            #'
            #' ## Dataset: `r name`
            #'
            #+ echo=FALSE
            finalTable <- colData

            #'
            #' ## Final sample table `r name`
            #'
            #+ savetable
            DT::datatable(finalTable, options=list(scrollX=TRUE))

            dim(finalTable)
            write_tsv(finalTable, file="!{group}.tsv")

            # run the version part
            cat(file="versions.yml", "!{task.process}:\naberrantexpression: 1.3.0")
        '''
}
