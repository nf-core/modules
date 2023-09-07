#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BTYPER3 } from '../../../../modules/nf-core/btyper3/main.nf'

workflow test_btyper3 {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BTYPER3 ( input )
}
