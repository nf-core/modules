#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_MASKFASTA } from '../../../../software/bedtools/maskfasta/main.nf' addParams( options: [:] )

workflow test_bedtools_maskfasta {
    bed   = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed", checkIfExists: true)
            ]
    fasta = [ file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true) ]

    BEDTOOLS_MASKFASTA ( bed, fasta )
}
