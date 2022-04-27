#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARDICTJAVA } from '../../../modules/vardictjava/main.nf'

workflow test_vardictjava {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    VARDICTJAVA ( input )
}
