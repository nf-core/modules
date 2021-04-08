#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUBBINS } from '../../../software/gubbins/main.nf' addParams( options: [:] )

workflow test_gubbins { 
    input = file(params.test_data['sarscov2']['genome']['alignment']['all_sites_alignment'], checkIfExists: true)

    GUBBINS ( input )
}
