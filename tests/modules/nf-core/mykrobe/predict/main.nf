#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MYKROBE_PREDICT } from '../../../../../modules/nf-core/mykrobe/predict/main.nf'

workflow test_mykrobe_predict {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    MYKROBE_PREDICT ( input, "typhi" )
}
