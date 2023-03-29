#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COREOGRAPH } from '../../../../modules/nf-core/coreograph/main.nf'

workflow test_coreograph {

    input = [
        [ id:'test', single_end:false ], // meta map
        file('/workspace/modules/tests/modules/nf-core/coreograph/TestdataCoreographHighRes.tif', checkIfExists: true)
    ]

    COREOGRAPH ( input )
}
