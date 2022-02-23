#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQKIT_PAIR } from '../../../../modules/seqkit/pair/main.nf'

workflow test_seqkit_pair {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) 
    ]

    SEQKIT_PAIR ( input )
}
