#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MUSCLE5_SUPER5 } from '../../../../../modules/nf-core/muscle5/super5/main.nf'

workflow test_muscle5_super5 {
    
    input = [
        [ id:'test' ], 
        fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    ]

    MUSCLE5_SUPER5 ( input )
}
