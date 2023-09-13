#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GEM3_GEM3INDEXER } from '../../../../../modules/nf-core/gem3/gem3indexer/main.nf'

workflow test_gem3_gem3indexer {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    GEM3_GEM3INDEXER ( input )
}
