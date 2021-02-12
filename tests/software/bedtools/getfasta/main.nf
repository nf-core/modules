#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_GETFASTA } from '../../../../software/bedtools/getfasta/main.nf' addParams( options: [:] )

workflow test_bedtools_getfasta {
    def bed   = [ file("${launchDir}/tests/data/bed/C.bed", checkIfExists: true) ]
    def fasta = [ file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true) ]

    BEDTOOLS_GETFASTA ( bed, fasta )
}
