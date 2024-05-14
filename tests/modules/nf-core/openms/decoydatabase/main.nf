#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { OPENMS_DECOYDATABASE } from '../../../../../modules/nf-core/openms/decoydatabase/main.nf'

workflow test_openms_decoydatabase {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['proteomics']['database']['yeast_ups'], checkIfExists: true)
    ]

    OPENMS_DECOYDATABASE ( input )
}
