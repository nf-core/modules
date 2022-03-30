#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KAIJU_KAIJU } from '../../../../modules/kaiju/kaiju/main.nf'

workflow test_kaiju_kaiju_single_end {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    db    = [
        file(params.test_data['sarscov2']['genome']['kaiju_fmi'], checkIfExists: true), // database
        file(params.test_data['sarscov2']['genome']['kaiju_nodes'], checkIfExists: true) // taxon nodes
    ]

    KAIJU_KAIJU ( input, db )
}

workflow test_kaiju_kaiju_paired_end {

    input = [
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
          file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
    ]
    db    = [
        file(params.test_data['sarscov2']['genome']['kaiju_fmi'], checkIfExists: true), // database
        file(params.test_data['sarscov2']['genome']['kaiju_nodes'], checkIfExists: true) // taxon nodes
    ]

    KAIJU_KAIJU ( input, db )
}
