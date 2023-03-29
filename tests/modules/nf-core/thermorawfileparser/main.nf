#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { THERMORAWFILEPARSER } from '../../../../modules/nf-core/thermorawfileparser/main.nf'

workflow test_thermorawfileparser {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['proteomics']['rawfile']['raw'], checkIfExists: true)
    ]

    THERMORAWFILEPARSER ( input )
}
