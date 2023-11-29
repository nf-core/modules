#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLTOOLS_INSULATION } from '../../../../../modules/nf-core/cooltools/insulation/main.nf'

workflow test_cooltools_insulation {


    input = [ [ id:'test' ], // meta map
              file("https://github.com/open2c/cooltools/raw/master/tests/data/CN.mm9.1000kb.cool", checkIfExists: true)
    ]

    COOLTOOLS_INSULATION ( input )
}
