#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_ZOOMIFY } from '../../../../modules/cooler/zoomify/main.nf' addParams( options: ['args':'-r 2,4,8', publish_files:[:]] )
include { COOLER_DUMP } from '../../../../modules/cooler/dump/main.nf' addParams( options: [:] )

workflow test_cooler_zoomify {
   input = [ [ id:'test' ], // meta map
            file(params.test_data['generic']['cooler']['test_merge_cool'], checkIfExists: true)]

    COOLER_ZOOMIFY ( input )
    COOLER_DUMP(COOLER_ZOOMIFY.out.mcool, "/resolutions/2")
}
