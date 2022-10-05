#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQOLD } from '../../../modules/fastqold/main.nf'

workflow test_fastqold {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    FASTQOLD ( input )
}
