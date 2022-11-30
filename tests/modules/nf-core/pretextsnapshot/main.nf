#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PRETEXTSNAPSHOT } from '../../../../modules/nf-core/pretextsnapshot/main.nf'

workflow test_pretextsnapshot {

    input = [
        [ id:'test' ], // meta map
        file("https://github.com/epaule/test-datasets/raw/pretextsnapshot/testdata/jaGalFasc40_2.pretext") // test data
    ]

    PRETEXTSNAPSHOT_ALL ( input )
}
