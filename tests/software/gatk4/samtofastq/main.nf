#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_SAMTOFASTQ } from '../../../../software/gatk4/samtofastq/main.nf' addParams( options: [:] )

workflow test_gatk4_samtofastq_single_end {

    def input = []
    input = [ [ id:'test', single_end: true ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/bam/test_single_end.bam", checkIfExists: true)] ]

    GATK4_SAMTOFASTQ ( input )
}

workflow test_gatk4_samtofastq_paired_end {

    def input = []
    input = [ [ id:'test', single_end: false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.bam", checkIfExists: true)] ]

    GATK4_SAMTOFASTQ ( input )
}
