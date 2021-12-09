#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SMOOTHXG } from '../../../modules/smoothxg/main.nf'

workflow test_smoothxg {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['assembly_gfa'], checkIfExists: true) ]
            ]

    SMOOTHXG ( input )
}
