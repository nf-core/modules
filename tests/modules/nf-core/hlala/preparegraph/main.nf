#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HLALA_PREPAREGRAPH } from '../../../../../modules/nf-core/hlala/preparegraph/main.nf'

workflow test_hlala_preparegraph {
    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['genome']['prg_input'], checkIfExists: true)
    ]
    HLALA_PREPAREGRAPH ( input )
}
