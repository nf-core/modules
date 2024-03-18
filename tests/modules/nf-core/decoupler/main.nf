#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DECOUPLER } from '../../../../modules/nf-core/decoupler/main.nf'

workflow test_decoupler {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    DECOUPLER ( input )
}
