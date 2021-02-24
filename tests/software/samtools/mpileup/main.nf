#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_MPILEUP } from '../../../../software/samtools/mpileup/main.nf' addParams( options: [:] )

workflow test_samtools_mpileup {

    def input = []
    def fasta = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test-sc2-artic-v3-sorted-trimmed.bam", checkIfExists: true) ]
    fasta = [  file("${launchDir}/tests/data/genomics/sarscov2/fasta/MN908947.3.fa", checkIfExists: true)  ]
    SAMTOOLS_MPILEUP ( input, fasta )
}
