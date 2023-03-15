process COUNTREADS {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each preprocess_data
        each bam_data

    output:
        tuple val(preprocess), val(groups), path("*.Rds")       , emit: result
        path "versions.yml"                                     , emit: versions

    shell:
        count_ranges = preprocess_data.countRanges
        preprocess = preprocess_data

        bam = bam_data.RNA_BAM_FILE
        groups = bam_data.DROP_GROUP

        countMode = bam_data.COUNT_MODE
        pairedEnd = bam_data.PAIRED_END
        strandSpec = bam_data.STRAND
        countOverlaps = bam_data.COUNT_OVERLAPS
        RNA_ID = bam_data.RNA_ID
        '''
            #!/usr/bin/env Rscript --vanilla

            suppressPackageStartupMessages({
                library(data.table)
                library(Rsamtools)
                library(BiocParallel)
                library(GenomicAlignments)
            })
            inputBAM <- "!{bam}"
            inputRanges <- "!{count_ranges}"

            # Get strand specific information from sample annotation
            sampleID <- "!{RNA_ID}"

            strand <- tolower("!{strandSpec}")
            count_mode <- "!{countMode}"
            paired_end <- as.logical("!{pairedEnd}")
            overlap <- as.logical("!{countOverlaps}")
            inter_feature <- ! overlap # inter_feature = FALSE does not allow overlapsi

            # infer preprocessing and strand info
            preprocess_reads <- NULL
            if (strand == "yes") {
                strand_spec <- T
            } else if (strand == "no"){
                strand_spec <- F
            } else if (strand == "reverse") {
                # set preprocess function for later
                preprocess_reads <- invertStrand
                strand_spec <- T
            } else {
                stop(paste("invalid strand information", strand))
            }

            # read files
            bam_file <- BamFile(inputBAM, yieldSize = 2e6)
            count_ranges <- readRDS(inputRanges)
            # set chromosome style
            seqlevelsStyle(count_ranges) <- seqlevelsStyle(bam_file)

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
                , BPPARAM = MulticoreParam(1)
            )
            colnames(se) <- sampleID
            saveRDS(se, paste0(sampleID, ".Rds"))

            # run the version part
            cat(file="versions.yml", "!{task.process}:\naberrantexpression: 1.3.0")
        '''
}
