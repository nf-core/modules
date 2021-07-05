#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ADAM_MARKDUPLICATES } from '../../../../software/adam/markduplicates/main.nf' addParams( options: [:] )

workflow test_adam_markduplicates {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ] ]

    ADAM_MARKDUPLICATES ( input )
}
