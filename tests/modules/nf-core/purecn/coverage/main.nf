#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURECN_COVERAGE } from '../../../../../modules/nf-core/purecn/coverage/main.nf'

workflow test_purecn_coverage {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    PURECN_COVERAGE ( input )
}
