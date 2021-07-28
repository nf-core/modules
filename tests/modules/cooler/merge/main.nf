#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_MERGE } from '../../../../modules/cooler/merge/main.nf' addParams( options: [:] )

workflow test_cooler_merge {

    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['homo_sapiens']['cooler']['test_merge_cool'], checkIfExists: true) ] ]

    COOLER_MERGE ( input )
}
