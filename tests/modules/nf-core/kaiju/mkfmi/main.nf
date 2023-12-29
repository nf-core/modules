#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KAIJU_MKFMI } from '../../../../../modules/nf-core/kaiju/mkfmi/main.nf'

workflow test_kaiju_mkfmi {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true)
    ]

    KAIJU_MKFMI ( input )
}
