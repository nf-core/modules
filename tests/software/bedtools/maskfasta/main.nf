#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_MASKFASTA } from '../../../../software/bedtools/maskfasta/main.nf' addParams( options: [:] )

workflow test_bedtools_maskfasta {
    def input,fasta = []
    bed = [ [ id:'test'],
              file("${launchDir}/tests/data/bed/C.bed", checkIfExists: true) ]
    fasta =  [ file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true) ]

    BEDTOOLS_MASKFASTA( bed, fasta )
}
