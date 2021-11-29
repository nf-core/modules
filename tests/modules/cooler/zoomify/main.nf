#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_ZOOMIFY } from '../../../../modules/cooler/zoomify/main.nf'

workflow test_cooler_zoomify {

   input = [
        [ id:'test' ], // meta map
        file(params.test_data['generic']['cooler']['test_merge_cool'], checkIfExists: true)
    ]

    COOLER_ZOOMIFY ( input )
}
