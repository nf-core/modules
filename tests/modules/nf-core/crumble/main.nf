#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CRUMBLE } from '../../../../modules/nf-core/crumble/main.nf'

workflow test_crumble {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    CRUMBLE ( input, [], [] )
}

workflow test_crumble_bedout {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    CRUMBLE ( input, [], true )
}
