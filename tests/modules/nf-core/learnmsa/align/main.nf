#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LEARNMSA_ALIGN } from '../../../../../modules/nf-core/learnmsa/align/main.nf'

workflow test_learnmsa_align {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['genome']['alignment']['informative_sites'], checkIfExists: true)
    ]

    LEARNMSA_ALIGN ( input )
}
