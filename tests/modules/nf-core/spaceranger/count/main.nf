#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPACERANGER_COUNT } from '../../../../../modules/nf-core/spaceranger/count/main.nf'

workflow test_spaceranger_count {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    SPACERANGER_COUNT ( input )
}
