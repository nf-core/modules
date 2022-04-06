#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR } from '../../../modules/untar/main.nf'

workflow test_untar {
    input = [
        [],
        file(params.test_data['sarscov2']['genome']['kraken2_tar_gz'], checkIfExists: true)
    ]

    UNTAR ( input )
}
