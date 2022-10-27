#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RTGTOOLS_PEDFILTER } from "$moduleDir/modules/nf-core/rtgtools/pedfilter/main.nf"

workflow test_rtgtools_pedfilter {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['justhusky_ped'], checkIfExists: true)
    ]

    RTGTOOLS_PEDFILTER ( input )
}
