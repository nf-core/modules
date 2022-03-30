#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { METAMAPS_MAPDIRECTLY } from '../../../../modules/metamaps/mapdirectly/main.nf'

workflow test_metamaps_mapdirectly {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['nanopore']['test2_fastq_gz'], checkIfExists: true)
    ]
    database = [
        file(params.test_data['sarscov2']['genome']['metamaps_db'], checkIfExists: true)
    ]

    METAMAPS_MAPDIRECTLY ( input )
}
