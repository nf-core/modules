#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_MULTIINTER } from '../../../../../modules/nf-core/bedtools/multiinter/main.nf'

workflow test_bedtools_multiinter {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
            file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)
        ]
    ]

    BEDTOOLS_MULTIINTER ( input, [] )
}

workflow test_bedtools_multiinter_genome {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
            file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)
        ]
    ]
    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)

    BEDTOOLS_MULTIINTER ( input, sizes )
}
