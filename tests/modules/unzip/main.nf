#!/usr/bin/env nextflow



include { UNZIP } from '../../../modules/unzip/main.nf'

workflow test_unzip {

    archive = [
        [],
        file(params.test_data['sarscov2']['genome']['ncbi_taxmap_zip'], checkIfExists: true)
    ]

    UNZIP ( archive )
}
