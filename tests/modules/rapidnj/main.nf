#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RAPIDNJ } from '../../../software/rapidnj/main.nf' addParams( options: [:] )

workflow test_rapidnj {
    
    input = [ file(params.test_data['sarscov2']['genome']['informative_sites_fas'], checkIfExists: true) ]

    RAPIDNJ ( input )
}
