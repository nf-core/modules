#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRANSDECODER_PREDICT } from '../../../../../modules/nf-core/transdecoder/predict/main.nf'
include { TRANSDECODER_LONGORF } from '../../../../../modules/nf-core/transdecoder/longorf/main.nf'

workflow test_transdecoder_predict {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)]
    
    TRANSDECODER_LONGORF ( input )
    TRANSDECODER_PREDICT ( input, TRANSDECODER_LONGORF.out.folder )

}

workflow test_transdecoder_predict_gz_input {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta_gz'], checkIfExists: true)]
    
    TRANSDECODER_LONGORF ( input )
    TRANSDECODER_PREDICT ( input, TRANSDECODER_LONGORF.out.folder )

}
