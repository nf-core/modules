#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MUSCLE } from '../../../modules/muscle/main.nf'
include { MUSCLE as MUSCLE_TREE } from '../../../modules/muscle/main.nf'

workflow test_muscle {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['genome']['informative_sites_fas'], checkIfExists: true) ]

    MUSCLE ( input )
    MUSCLE_TREE ( MUSCLE.out.aligned_fasta )
}
