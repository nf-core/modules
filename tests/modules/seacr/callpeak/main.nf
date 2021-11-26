#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEACR_CALLPEAK } from '../../../../modules/seacr/callpeak/main.nf'

workflow test_seacr_callpeak {
    input = [ [ id:'test_1'],
              file(params.test_data['homo_sapiens']['illumina']['cutandrun_bedgraph_test_1'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['cutandrun_bedgraph_test_2'], checkIfExists: true) 
            ]

    SEACR_CALLPEAK ( input, 0.05 )
}

workflow test_seacr_callpeak_threshold {
    input = [ [ id:'test_1'],
              file(params.test_data['homo_sapiens']['illumina']['cutandrun_bedgraph_test_1'], checkIfExists: true),
              []
            ]

    SEACR_CALLPEAK ( input, 0.05 )
}