#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GFAFFIX } from '../../../../modules/nf-core/gfaffix/main.nf'

workflow test_gfaffix {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['assembly_gfa'], checkIfExists: true)]
            ]

    GFAFFIX ( input )
}
