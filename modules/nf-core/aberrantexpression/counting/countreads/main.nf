process COUNTREADS {
//    tag "$meta.id"
    label 'process_medium'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        tuple path(txdb), path(gene_map), path(ranges)
        each path(bam)

    output:
        path("*.Rds"), emit: counts
        path "versions.yml"           , emit: versions
    
    shell:
    '''
    #!/usr/bin/env Rscript --vanilla
    suppressPackageStartupMessages({
      library(data.table)
      library(Rsamtools)
      library(BiocParallel)
      library(GenomicAlignments)
    })
    inputBAM <- "!{bam}"
    inputRanges <- "!{ranges}"

    # Get strand specific information from sample annotation
    sampleID <- gsub(".bam", "", inputBAM)
    
    # TODO parse the input/parameters 
    # check: https://github.com/gagneurlab/drop/blob/37a91d34ff0d8ef51906ae7edb004cdf53c846f0/drop/modules/aberrant-expression-pipeline/Counting/countReads.R#L31
    
    # read files
    bam_file <- BamFile(inputBAM, yieldSize = 1e5)
    count_ranges <- readRDS(inputRanges)
    
    # set chromosome style
    seqlevelsStyle(count_ranges) <- seqlevelsStyle(bam_file)
    
    count_mode <- "Union"
    paired_end <- TRUE
    strand_spec <- TRUE
    overlap <- FALSE
    inter_feature <- ! overlap # inter_feature = FALSE does not allow overlapsi
    preprocess_reads <- NULL

    # start counting
    se <- summarizeOverlaps(
          count_ranges
        , bam_file
        , mode = count_mode
        , singleEnd = !paired_end
        , ignore.strand = !strand_spec  # FALSE if done strand specifically
        , fragments = F
        , count.mapped.reads = T
        , inter.feature = inter_feature # TRUE: reads mapping to multiple features are dropped
        , preprocess.reads = preprocess_reads
    )
    colnames(se) <- sampleID
    saveRDS(se, paste0(sampleID, ".Rds"))
    
    # run the version part
    cat(file="versions.yml", "!{task.process}:\naberrantexpression: 1.3.0")
    '''
}
