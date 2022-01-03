#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BANDAGE_IMAGE } from '../../../../modules/bandage/image/main.nf'

workflow test_bandage_image {
    input = [
        [ id:'B-3106' ], // meta map
        file( params.test_data['sarscov2']['illumina']['assembly_gfa'], checkIfExists: true)
    ]

    BANDAGE_IMAGE ( input )
}
