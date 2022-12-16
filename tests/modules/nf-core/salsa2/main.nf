#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SALSA2 } from '../../../../modules/nf-core/salsa2/main.nf'

workflow test_salsa2 {

    input = [
        [ id:'test', single_end: false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]
    bed = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)

    SALSA2 ( input, bed, [], [], [] )
}
