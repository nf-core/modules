#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KAIJU } from '../../../modules/kaiju/main.nf'

workflow test_kaiju_single_end {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    db    = [
        file(params.test_data['sarscov2']['genome']['kaiju_fmi'], checkIfExists: true), // database
        file(params.test_data['sarscov2']['genome']['kaiju_nodes'], checkIfExists: true), // taxon nodes
        file(params.test_data['sarscov2']['genome']['kaiju_names'], checkIfExists: true) // taxon names
    ]
    taxon_rank = "species"

    KAIJU ( input, db, taxon_rank )
}

workflow test_kaiju_paired_end {

    input = [
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
          file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
    ]
    db    = [
        file(params.test_data['sarscov2']['genome']['kaiju_fmi'], checkIfExists: true), // database
        file(params.test_data['sarscov2']['genome']['kaiju_nodes'], checkIfExists: true), // taxon nodes
        file(params.test_data['sarscov2']['genome']['kaiju_names'], checkIfExists: true) // taxon names
    ]
    taxon_rank = "species"

    KAIJU ( input, db, taxon_rank )
}
