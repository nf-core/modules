#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CDHIT } from '../../../../modules/nf-core/cdhit/main.nf'

workflow test_cdhit {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    CDHIT ( input )
}
