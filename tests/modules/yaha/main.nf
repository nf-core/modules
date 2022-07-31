#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { YAHA } from '../../../modules/yaha/main.nf'

workflow test_yaha {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    YAHA ( input )
}
