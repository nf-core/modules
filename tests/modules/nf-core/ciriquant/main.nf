#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CIRIQUANT } from '../../../../modules/nf-core/ciriquant/main.nf'

workflow test_ciriquant {

    reads = [
        [ id:'fust1_1', single_end:false ],
        file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/fastq/fust1_rep1_1.fastq.gz"),
        file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/fastq/fust1_rep1_2.fastq.gz")
    ]

    gtf = [ file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/reference/chrI.gtf") ]

    fasta = [ file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/reference/chrI.fa") ]

    bwa = [ file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/bwa") ]

    hisat2 = [ file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/hisat2") ]


    CIRIQUANT ( reads, gtf, fasta, bwa, hisat2 )
}
