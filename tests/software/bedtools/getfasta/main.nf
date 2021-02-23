#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_GETFASTA } from '../../../../software/bedtools/getfasta/main.nf' addParams( options: [:] )

workflow test_bedtools_getfasta {
    def bed   = [ file("${launchDir}/tests/data/genomics/sarscov2/bed/sarscov2.bed", checkIfExists: true) ]
    def fasta = [ file("${launchDir}/tests/data/genomics/sarscov2/fasta/GCA_011545545.1_ASM1154554v1_genomic.fna", checkIfExists: true) ]

    BEDTOOLS_GETFASTA ( bed, fasta )
}
