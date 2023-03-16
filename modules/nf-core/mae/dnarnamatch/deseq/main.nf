process DESEQ_QC {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each allelic_counts

    output:
        tuple val(vcf), val(rnaID), path("output.Rds")     , emit: result

    shell:
        vcf     = allelic_counts.vcf
        rnaID   = allelic_counts.rnaID

        counts  = allelic_counts.counts
        '''
            #!/usr/bin/env Rscript --vanilla

            suppressPackageStartupMessages({
                library(GenomicRanges)
                library(tMAE)
            })

            # Read MA counts for qc
            qc_counts <- fread("!{counts}", fill=TRUE)
            # qc_counts <- qc_counts[!is.na(position)]

            # Run DESeq
            rmae <- DESeq4MAE(qc_counts, minCoverage = 10)
            rmae[, RNA_GT := '0/1']
            rmae[altRatio < .2, RNA_GT := '0/0']
            rmae[altRatio > .8, RNA_GT := '1/1']
            rmae[, position := as.numeric(position)]

            # Convert to granges
            qc_gr <- GRanges(seqnames = rmae$contig,
                            ranges = IRanges(start = rmae$position, end = rmae$position),
                            strand = '*')
            mcols(qc_gr) = DataFrame(RNA_GT = rmae$RNA_GT)

            saveRDS(qc_gr, "output.Rds")
        '''
}
