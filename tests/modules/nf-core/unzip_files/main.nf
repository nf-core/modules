#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNZIP_FILES } from '../../../../modules/nf-core/unzip_files/main.nf'

workflow test_unzip_files {

    archive = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['ncbi_taxmap_zip'], checkIfExists: true)
    ]

    UNZIP_FILES ( archive )
}

