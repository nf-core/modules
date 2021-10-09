#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SRATOOLS_FASTERQDUMP } from '../../../../modules/sratools/fasterqdump/main.nf' addParams( options: [:] )

workflow test_sratools_fasterqdump {

    input = [
        [ id:'test', single_end:false ], // meta map
        file('tests/modules/sratools/fasterqdump/ERR2815334', checkIfExists: true)
    ]

    SRATOOLS_FASTERQDUMP ( input )
}
