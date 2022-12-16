#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CDHIT } from '../../../../modules/nf-core/cdhit/cdhit/main.nf'

workflow test_cdhit {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['proteomics']['database']['yeast_UPS_fasta'], checkIfExists: true)
    ]

    CDHIT ( input )
}
