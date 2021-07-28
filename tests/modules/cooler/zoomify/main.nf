#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_ZOOMIFY } from '../../../../modules/cooler/zoomify/main.nf' addParams( options: ['args':'-r 2,4,8'] )

workflow test_cooler_zoomify {
   input = [ [ id:'test' ], // meta map
            file(params.test_data['homo_sapiens']['cooler']['test_merge_cool'], checkIfExists: true)]

    COOLER_ZOOMIFY ( input )
}
