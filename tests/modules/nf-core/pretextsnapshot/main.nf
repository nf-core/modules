#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PRETEXTSNAPSHOT as PRETEXTSNAPSHOT_ALL} from '../../../../modules/nf-core/pretextsnapshot/main.nf'

workflow test_pretextsnapshot {

    input = [
        [ id:'test' ], // meta map
        file("/tmp/jaGalFasc40_2.pretext") // test data (copy it locally until in test-data)
    ]

    PRETEXTSNAPSHOT_ALL ( input )
}
