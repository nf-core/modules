#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRANSDECODER_PREDICT } from '../../../../modules/transdecoder/predict/main.nf'

workflow test_transdecoder_predict {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    ]

    TRANSDECODER_PREDICT ( input )
}
