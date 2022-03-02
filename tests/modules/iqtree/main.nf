#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IQTREE } from '../../../modules/iqtree/main.nf'

workflow test_iqtree {
    
    input = [ file(params.test_data['sarscov2']['genome']['informative_sites_fas'], checkIfExists: true) ]

    IQTREE ( input, '' )
}
