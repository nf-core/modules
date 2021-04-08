#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNPSITES } from '../../../software/snpsites/main.nf' addParams( options: [:] )

workflow test_snpsites {
    
    input = file(params.test_data['sarscov2']['genome']['alignment']['all_sites_alignment'], checkIfExists: true)

    SNPSITES ( input )
}
