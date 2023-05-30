#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IPHOP_PREDICT } from '../../../../../modules/nf-core/iphop/predict/main.nf'

workflow test_iphop_predict {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    IPHOP_PREDICT ( input )
}
