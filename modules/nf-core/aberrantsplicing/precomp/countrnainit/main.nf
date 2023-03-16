process COUNTRNA_INIT {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each groupData

    output:
        tuple val(group), \
            path("FRASER_output/savedObjects/raw-local-*/fds-object.RDS"), \
            path("FRASER_output")   , emit: result

    shell:
        group = groupData.group
        tsv   = groupData.tsv
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

            dataset     <- "!{group}"
            colDataFile <- "!{tsv}"

            # Create initial FRASER object
            col_data <- fread(colDataFile)

            fds <- FraserDataSet(colData = col_data,
                                name       = paste0("raw-local-", dataset))

            # Add paired end and strand specificity to the fds
            pairedEnd(fds) <- colData(fds)$PAIRED_END
            strandSpecific(fds) <- 'no'
            if(uniqueN(colData(fds)$STRAND) == 1){
            strandSpecific(fds) <- unique(colData(fds)$STRAND)
            }

            # Save initial FRASER dataset
            fds <- saveFraserDataSet(fds)

            message(date(), ": FRASER object initialized for ", dataset)

            # run the version part
            cat(file="versions.yml", "!{task.process}:\naberrantexpression: 1.3.0")
        '''
}
