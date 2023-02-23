#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURGEDUPS_SPLITFA } from '../../../../../modules/nf-core/purgedups/splitfa/main.nf'

workflow test_purgedups_splitfa {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    PURGEDUPS_SPLITFA ( input )
}
