#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { OPENMS_IDMERGER } from '../../../../../modules/nf-core/openms/idmerger/main.nf'

workflow test_openms_idmerger {

    input = [
        [ id:'test' ], // meta map
        [
            file(params.test_data['proteomics']['database']['idxml1'], checkIfExists: true),
            file(params.test_data['proteomics']['database']['idxml2'], checkIfExists: true)
        ]
    ]

    OPENMS_IDMERGER ( input )
}
