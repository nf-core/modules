#!/usr/bin/env nextflow



include { BAKTA } from '../../../modules/bakta/main.nf'

workflow test_bakta {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BAKTA ( input, [], [], [] )
}
