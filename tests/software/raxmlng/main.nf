#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RAXMLNG as RAXMLNG_NO_BOOTSTRAP } from '../../../software/raxmlng/main.nf' addParams( options: [args:'--model GTR+G']                 )
include { RAXMLNG as RAXMLNG_BOOTSTRAP    } from '../../../software/raxmlng/main.nf' addParams( options: [args:'--model GTR+G --bs-trees 1000'] )

/*
 * Test without bootstrapping
 */

workflow test_raxmlng_no_bootstrap {
    
    input = [ file(params.test_data['sarscov2']['genome']['informative_sites_fas'], checkIfExists: true) ]

    RAXMLNG_NO_BOOTSTRAP ( input )
}

/*
 * Test with bootstrapping
 */

workflow test_raxmlng_bootstrap {
    
    input = [ file(params.test_data['sarscov2']['genome']['informative_sites_fas'], checkIfExists: true) ]

    RAXMLNG_BOOTSTRAP ( input )
}
