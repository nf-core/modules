#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MYKROBE_PREDICT } from '../../../../modules/mykrobe/predict/main.nf'
include { GENOTYPHI_PARSE } from '../../../../modules/genotyphi/parse/main.nf'

workflow test_mykrobe_predict {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    MYKROBE_PREDICT ( input, "typhi" )
    GENOTYPHI_PARSE ( MYKROBE_PREDICT.out.json )
}
