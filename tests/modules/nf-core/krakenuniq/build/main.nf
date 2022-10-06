#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRAKENUNIQ_BUILD } from '../../../../../modules/nf-core/krakenuniq/build/main.nf'

workflow test_krakenuniq_build {


    input = [
        <lib_dir>,
        <tax_dir>,
        file(params.test_data['sarscov2']['metagenome']['krona_taxonomy'], checkIfExists: true)
    ]

    KRAKENUNIQ_BUILD ( input )
}
