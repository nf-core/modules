#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_MPILEUP } from '../../../../software/bcftools/mpileup/main.nf' addParams( options: ['args2': '--no-version --ploidy 1 --multiallelic-caller',
                                                                                                       'args3': '--no-version' ] )

workflow test_bcftools_mpileup {
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true) ]]
    fasta = [ file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true) ]

    BCFTOOLS_MPILEUP ( input, fasta )
}
