#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_MPILEUP } from '../../../../software/samtools/mpileup/main.nf' addParams( options: [:] )

workflow test_samtools_mpileup {

    def input = []
    def fasta = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true) ]
    fasta = [  file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)  ]
    SAMTOOLS_MPILEUP ( input, fasta )
}
