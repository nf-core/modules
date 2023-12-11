#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { OPENMS_IDFILTER } from '../../../../../modules/nf-core/openms/idfilter/main.nf'

workflow test_openms_idfilter {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['proteomics']['openms']['idxml1'], checkIfExists: true),
        file(params.test_data['proteomics']['openms']['idxml2'], checkIfExists: true)
    ]

    OPENMS_IDFILTER ( input )
}
