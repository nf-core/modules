#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARACLU } from '../../../modules/paraclu/main.nf'

workflow test_paraclu {

    input = [[  id:'test' ], // meta map
                file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
                ]
    min_cluster = 30

    PARACLU ( input, min_cluster )
}
