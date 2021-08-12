#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GSTAMA_MERGE } from '../../../../modules/gstama/merge/main' addParams( options: [:] )

workflow test_gstama_merge {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['galgal6']['bed']['chr30'], checkIfExists: true),
            file(params.test_data['galgal6']['bed']['chr31'], checkIfExists: true)
        ]
    ]

    GSTAMA_MERGE ( input )
}
