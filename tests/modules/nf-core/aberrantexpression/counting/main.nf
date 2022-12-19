#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREPROCESSGENEANNOTATION  } from '../../../../../modules/nf-core/aberrantexpression/counting/preprocessgeneannotation'
include { COUNTREADS                } from '../../../../../modules/nf-core/aberrantexpression/counting/countreads'

process prepare_data {
    output:
        path "drop_demo_data-main/Data/rna_bam/*bam"
        file "drop_demo_data-main/Data/*gtf"
    
    script:
    """
        curl -L https://github.com/nickhsmith/drop_demo_data/archive/refs/heads/main.zip -o data.zip && \
            unzip data
    """
}

process merge_data {
    conda (params.enable_conda ? "TODO" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"
    input:
        path(counts)
    output:
        file("merged.counts.Rds")
    shell:
    """
    #!/usr/bin/env Rscript --vanilla
    library(SummarizedExperiment)
    library(GenomicRanges)
    
    inputFiles <- "!{counts}"
    
    # TODO: adapt the code from here: https://github.com/gagneurlab/drop/blob/37a91d34ff0d8ef51906ae7edb004cdf53c846f0/drop/modules/aberrant-expression-pipeline/Counting/mergeCounts.R#L20
    files <- strsplit(inputFiles, " ")[[1]]
    rds <- lapply(files, readRDS)
    colD <- do.call(rbind, lapply(rds, colData))
    assC <- do.call(cbind, lapply(rds, assay, "counts"))
    
    se <- SummarizedExperiment(
        rowRanges=rowRanges(rds[[1]]),
        assays=SimpleList(counts=assC),
        colData=colD)
    
    # save in RDS format
    saveRDS(se, "merged.counts.Rds")
    """
}

workflow test_aberrantexpression_counting {
    // get all required input data from demo into corresponding channels  
    (ch_bam, ch_gtf) = prepare_data()
    
    // run the preprocessing only on the gtf channel
    PREPROCESSGENEANNOTATION(ch_gtf)
    
    // run the counting per BAM file
    COUNTREADS(PREPROCESSGENEANNOTATION.out.count_ranges, ch_bam.flatten() )

    // merge all BAM files into a single one
    merge_data( COUNTREADS.out.counts.collect() )
}
