#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MUSCLE5_SUPER5 } from '../../../../../modules/nf-core/muscle5/super5/main.nf'

workflow test_muscle5_super5 {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    MUSCLE5_SUPER5 ( input )
}
