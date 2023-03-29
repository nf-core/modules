#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNZIPFILES } from '../../../../modules/nf-core/unzipfiles/main.nf'

workflow test_unzipfiles {

    archive = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['ncbi_taxmap_zip'], checkIfExists: true)
    ]

    UNZIPFILES ( archive )
}

