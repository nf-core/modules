#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GAWK } from '../../../../modules/nf-core/gawk/main.nf'

workflow test_gawk {

    // A small example of a conversion of fasta index to bed format

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    GAWK ( input )
}
