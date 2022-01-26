#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SVDB_MERGE } from '../../../../modules/svdb/merge/main.nf'

workflow test_svdb_merge {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) 
    ]

    SVDB_MERGE ( input )
}
