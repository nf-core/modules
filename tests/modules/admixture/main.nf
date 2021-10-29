#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ADMIXTURE } from '../../../modules/admixture/main.nf' addParams( options: [:] )

workflow test_admixture {
    
    input = [ [ id:'test', single_end:false ], // meta map
             file(params.test_data['homosapiens']['genome']['genome_bed'], checkIfExists: true) ]
    ADMIXTURE ( input )
}
