process DEFINE_DATASET_FROM_ANNO {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        path(sampleAnnotation)
        path(bams)

    output:
        path("fraser.tsv")            , emit: fraser

    shell:
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
    })

    #+ input
    outFile       <- "fraser.tsv"
    annoFile      <- "!{sampleAnnotation}"
    fileMapFile   <- "!{bams}"

    #+ dataset name

    anno    <- fread(annoFile)
    mapping <- strsplit(fileMapFile, " ") %>%
        map(function(x) c(str_split(x, ".")[[1]], x))
    print(mapping)

    setnames(annoSub, "RNA_ID", "sampleID")
    annoSub <- anno[RNA_ID]
    colData <- merge(annoSub,
        mapping[FILE_TYPE == "RNA_BAM_FILE", .(sampleID=ID, bamFile=FILE_PATH)])
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
    write_tsv(finalTable, file=outFile)

    # run the version part
    cat(file="versions.yml", "!{task.process}:\naberrantexpression: 1.3.0")
    '''
}
