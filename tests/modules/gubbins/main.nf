#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUBBINS } from '../../../modules/gubbins/main.nf'

workflow test_gubbins { 
    input = file(params.test_data['sarscov2']['genome']['all_sites_fas'], checkIfExists: true)

    GUBBINS ( input )
}
