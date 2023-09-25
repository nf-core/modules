#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GEM3_GEM3MAPPER } from '../../../../../modules/nf-core/gem3/gem3mapper/main.nf'

workflow test_gem3_gem3mapper {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    GEM3_GEM3MAPPER ( input )
}
