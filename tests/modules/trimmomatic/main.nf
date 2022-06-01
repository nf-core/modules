#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRIMMOMATIC } from '../../../modules/trimmomatic/main.nf'

workflow test_trimmomatic {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    TRIMMOMATIC ( input )
}
