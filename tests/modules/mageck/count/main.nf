#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAGECK_COUNT } from '../../../../modules/mageck/count/main.nf'

workflow test_mageck_count {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    MAGECK_COUNT ( input )
}
