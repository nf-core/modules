#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LOFREQ3_PREPROCESSING } from '../../../../../modules/nf-core/lofreq3/preprocessing/main.nf'

workflow test_lofreq3 {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    LOFREQ3 ( input )
}
