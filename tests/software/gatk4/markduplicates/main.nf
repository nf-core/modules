#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MARKDUPLICATES } from '../../../../software/gatk4/markduplicates/main.nf' addParams( options: [:] )

workflow test_gatk4_markduplicates {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true) ]

    GATK4_MARKDUPLICATES ( input )
}
