#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARACLU } from '../../../modules/paraclu/main.nf' addParams( options: [args:30] )

workflow test_paraclu {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) ]

    PARACLU ( input )
}
