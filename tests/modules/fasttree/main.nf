#!/usr/bin/env nextflow



include { FASTTREE } from '../../../modules/fasttree/main.nf'

workflow test_fasttree {
    
    input = [ file(params.test_data['sarscov2']['genome']['informative_sites_fas'], checkIfExists: true) ]

    FASTTREE ( input )
}
