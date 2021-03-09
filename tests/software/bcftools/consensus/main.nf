#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_CONSENSUS } from '../../../../software/bcftools/consensus/main.nf' addParams( options: [:] )

workflow test_bcftools_consensus {
    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/vcf/test.vcf.gz", checkIfExists: true) ],
              [ file("${launchDir}/tests/data/genomics/sarscov2/vcf/test.vcf.gz.tbi", checkIfExists: true) ],
              [ file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) ] ]

    BCFTOOLS_CONSENSUS ( input )
}
