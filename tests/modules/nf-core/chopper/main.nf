#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP } from '../../../../modules/nf-core/gunzip/main.nf'
include { CHOPPER } from '../../../../modules/nf-core/chopper/main.nf'

workflow test_chopper {

    input = [
        [id:'test_out' ], // meta map
        file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true)
    ]

    CHOPPER ( input )
}
