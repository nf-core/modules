#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_BGZIP as BCFTOOLS_BGZIP } from '../../../../software/bcftools/bgzip/main.nf' addParams( options: [:] )
include { BCFTOOLS_TABIX as BCFTOOLS_TABIX } from '../../../../software/bcftools/tabix/main.nf' addParams( options: [:] )
include { BCFTOOLS_CONSENSUS as BCFTOOLS_CONSENSUS } from '../../../../software/bcftools/consensus/main.nf' addParams( options: [:] )

workflow test_bcftools_consensus {

    def input = []
    input = [ [ id:'test' ], // meta map
            file("${launchDir}/tests/data/vcf/test.vcf", checkIfExists: true) ]
    def fasta = []
    fasta = Channel.of ([ [ id:'test' ], // meta map
            file("${launchDir}/tests/data/vcf/test.consensus.fa", checkIfExists: true) ])

    BCFTOOLS_BGZIP ( input )
    BCFTOOLS_TABIX ( BCFTOOLS_BGZIP.out.vcf )
    BCFTOOLS_CONSENSUS ( BCFTOOLS_BGZIP.out.vcf
                                    .join( BCFTOOLS_TABIX.out.tbi, by: [0] )
                                    .join( fasta, by: [0] ) )
}
