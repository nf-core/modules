#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GPROFILER2_GOST } from '../../../../../modules/nf-core/gprofiler2/gost/main.nf'

workflow test_gprofiler2_gost {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['mus_musculus']['genome']['deseq_results'], checkIfExists: true)
    ]

    GPROFILER2_GOST ( input )
}
