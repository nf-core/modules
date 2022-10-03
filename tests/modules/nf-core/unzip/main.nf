#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNZIP } from '../../../modules/unzip/main.nf'

workflow test_unzip {

    archive = [
        [],
        file(params.test_data['sarscov2']['genome']['ncbi_taxmap_zip'], checkIfExists: true)
    ]

    UNZIP ( archive )
}
