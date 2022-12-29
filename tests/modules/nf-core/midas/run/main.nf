#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MIDAS_RUN } from '../../../../../modules/nf-core/midas/run/main.nf'

workflow test_midas_run {
    
    input = [ 
        [ id:'test', single_end:true ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
    ]

    MIDAS_RUN ( input, [], "species")
}
