#!/usr/bin/env nextflow



include { SNPSITES } from '../../../modules/snpsites/main.nf'

workflow test_snpsites {
    
    input = file(params.test_data['sarscov2']['genome']['all_sites_fas'], checkIfExists: true)

    SNPSITES ( input )
}
