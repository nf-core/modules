#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_MPILEUP } from '../../../../software/bcftools/mpileup/main.nf' addParams( options: ['args2': '--no-version --ploidy 1 --multiallelic-caller'] )

workflow test_bcftools_mpileup {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/bam/test-sc2-artic-v3.bam", checkIfExists: true) ]]
    fasta = [ file("${launchDir}/tests/data/fasta/sarscov2/MN908947.3.fa", checkIfExists: true) ]

    BCFTOOLS_MPILEUP ( input, fasta )

}
