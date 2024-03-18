#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_ZOOMIFY } from '../../../../../modules/nf-core/cooler/zoomify/main.nf'
include { COOLER_DUMP    } from '../../../../../modules/nf-core/cooler/dump/main.nf'

workflow test_cooler_zoomify {

   input = [
        [ id:'test' ], // meta map
        file(params.test_data['generic']['cooler']['test_merge_cool'], checkIfExists: true)
    ]

    COOLER_ZOOMIFY ( input )
    COOLER_DUMP ( COOLER_ZOOMIFY.out.mcool.combine([2,4,8]) )
}
