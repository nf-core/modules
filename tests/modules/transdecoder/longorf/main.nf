#!/usr/bin/env nextflow



include { TRANSDECODER_LONGORF } from '../../../../modules/transdecoder/longorf/main.nf'

workflow test_transdecoder_longorf {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) 
    ]

    TRANSDECODER_LONGORF ( input )
}
