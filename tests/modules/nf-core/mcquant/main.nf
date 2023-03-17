#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MCQUANT } from '../../../../modules/nf-core/mcquant/main.nf'

workflow test_mcquant {

    image = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    mask = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    markerfile = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    MCQUANT ( image, mask, markerfile )
}
