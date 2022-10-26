#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARIBA_GETREF } from '../../../../../modules/nf-core/ariba/getref/main.nf'
include { ARIBA_RUN } from '../../../../../modules/nf-core/ariba/run/main.nf'

workflow test_ariba_run {
    
    input = [ [ id:'test', single_end:false ], // meta map
              [  file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                 file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    ARIBA_GETREF ( "card" )
    ARIBA_RUN ( input, ARIBA_GETREF.out.db)
}
