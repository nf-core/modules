#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RAXMLNG } from '../../../software/raxmlng/main.nf' addParams( options: [:] )

/*
 * Test without bootstrapping
 */

workflow test_raxmlng_no_bootstrap {
    
    input = [ file(params.test_data['sarscov2']['genome']['informative_sites_fas'], checkIfExists: true) ]

    RAXMLNG ( input )
}

/*
 * Test with bootstrapping
 */

workflow test_raxmlng_bootstrap {
    
    input = [ file(params.test_data['sarscov2']['genome']['informative_sites_fas'], checkIfExists: true) ]

    RAXMLNG ( input )
}
