#!/usr/bin/env nextflow



include { UNTAR } from '../../../modules/untar/main.nf'

workflow test_untar {
    input = [
        [],
        file(params.test_data['sarscov2']['genome']['kraken2_tar_gz'], checkIfExists: true)
    ]

    UNTAR ( input )
}
