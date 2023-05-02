#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NEXTDENOVO } from '../../../../modules/nf-core/nextdenovo/main.nf'

workflow test_nextdenovo {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true)
    ]

    NEXTDENOVO ( input, "hifi", "100k" )
}
