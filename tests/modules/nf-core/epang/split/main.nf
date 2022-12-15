#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EPANG_SPLIT } from '../../../../../modules/nf-core/epang/split/main.nf'

workflow test_epang_split {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    EPANG_SPLIT ( input )
}
