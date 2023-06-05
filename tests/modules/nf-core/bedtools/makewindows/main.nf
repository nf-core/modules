#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_MAKEWINDOWS } from '../../../../../modules/nf-core/bedtools/makewindows/main.nf'

workflow test_bedtools_makewindows_bed {

    input = [
        [ id:'test2'],
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    ]

    BEDTOOLS_MAKEWINDOWS ( input )
}

workflow test_bedtools_makewindows_fai {

    input = [
        [ id:'test2'],
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    BEDTOOLS_MAKEWINDOWS ( input )
}
