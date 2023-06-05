#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { WGSIM } from '../../../../modules/nf-core/wgsim/main.nf'

workflow test_wgsim {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    WGSIM ( input )
}
