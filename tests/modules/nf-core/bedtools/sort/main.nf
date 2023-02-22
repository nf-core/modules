#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_SORT } from '../../../../../modules/nf-core/bedtools/sort/main.nf'

workflow test_bedtools_sort {
    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
            ]

    BEDTOOLS_SORT ( input, [] )
}

workflow test_bedtools_sort_with_genome {
    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
            ]
    fai = Channel.of(file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true))

    BEDTOOLS_SORT ( input, [] )
}