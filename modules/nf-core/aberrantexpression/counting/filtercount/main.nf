process FILTERCOUNT {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each mergecount_data
        val params

    output:
        tuple val(preprocess), val(groupName), path("ods.Rds")      , emit: result
        path "versions.yml"                                         , emit: versions

    shell:
        txtDB = mergecount_data.preprocess.txtDb
        total_count = mergecount_data.totalCounts

        groupName = mergecount_data.groupName

        fpkmCutoff = params.fpkmCutoff

        preprocess = mergecount_data.preprocess
        '''
            #!/usr/bin/env Rscript --vanilla

            suppressPackageStartupMessages({
                library(data.table)
                library(GenomicFeatures)
                library(SummarizedExperiment)
                library(OUTRIDER)
            })

            counts <- readRDS("!{total_count}")
            ods <- OutriderDataSet(counts)
            txdb <- loadDb("!{txtDB}")

            # filter not expressed genes
            fpkmCutoff <- !{fpkmCutoff}
            ods <- filterExpression(ods, gtfFile=txdb, filter=FALSE,
                                    fpkmCutoff, addExpressedGenes=TRUE)

            # add column for genes with at least 1 gene
            rowData(ods)$counted1sample = rowSums(assay(ods)) > 0

            # External data check
            if (is.null(ods@colData$GENE_COUNTS_FILE)){ #column does not exist in sample annotation table
                has_external <- FALSE
            }else if(all(is.na(ods@colData$GENE_COUNTS_FILE))){ #column exists but it has no values
                has_external <- FALSE
            }else if(all(ods@colData$GENE_COUNTS_FILE == "")){ #column exists with non-NA values but this group has all empty strings
                has_external <- FALSE
            }else{ #column exists with non-NA values and this group has at least 1 non-empty string
                has_external <- TRUE
            }

            if(has_external){
                ods@colData$isExternal <- as.factor(ods@colData$GENE_COUNTS_FILE != "")
            }else{
                ods@colData$isExternal <- as.factor(FALSE)
            }


            # Save the ods before filtering to preserve the original number of genes
            saveRDS(ods, "ods.Rds")

            # run the version part
            cat(file="versions.yml", "!{task.process}:\naberrantexpression: 1.3.0")
        '''
}
