#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAPAD_INDEX } from '../../../../../modules/nf-core/mapad/index/main.nf'

workflow test_mapad_index {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    MAPAD_INDEX ( input )
}
