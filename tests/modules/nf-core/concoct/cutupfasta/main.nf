#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CONCOCT_CUTUPFASTA } from '../../../../../modules/nf-core/concoct/cutupfasta/main.nf'

workflow test_concoct_cutupfasta {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    produce_bedfile = true

    CONCOCT_CUTUPFASTA ( input, produce_bedfile )
}
