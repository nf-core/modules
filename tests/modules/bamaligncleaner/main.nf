#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAMALIGNCLEANER } from '../../../modules/bamaligncleaner/main.nf'

workflow test_bamaligncleaner {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true) ]

    BAMALIGNCLEANER ( input )
}
