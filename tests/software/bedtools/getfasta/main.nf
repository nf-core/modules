#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_GETFASTA } from '../../../../software/bedtools/getfasta/main.nf' addParams( options: [:] )

workflow test_bedtools_getfasta {
    def bed   = [ file("${launchDir}/tests/data/genomics/sarscov2/bed/test.bed", checkIfExists: true) ]
    def fasta = [ file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) ]

    BEDTOOLS_GETFASTA ( bed, fasta )
}
