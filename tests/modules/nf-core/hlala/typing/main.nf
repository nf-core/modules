#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HLALA_TYPING } from '../../../../../modules/nf-core/hlala/typing/main.nf'

workflow test_hlala_typing {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    HLALA_TYPING ( input )
}
