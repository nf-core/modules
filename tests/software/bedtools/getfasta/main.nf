#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_GETFASTA } from '../../../../software/bedtools/getfasta/main.nf' addParams( options: [:] )

workflow test_bedtools_getfasta {
    bed   = [ file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed", checkIfExists: true) ]
    fasta = [ file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true) ]

    BEDTOOLS_GETFASTA ( bed, fasta )
}
