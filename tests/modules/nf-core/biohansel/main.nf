#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BIOHANSEL as BIOHANSEL_PE } from '../../../../modules/nf-core/biohansel/main.nf'
include { BIOHANSEL as BIOHANSEL_SE } from '../../../../modules/nf-core/biohansel/main.nf'
include { BIOHANSEL as BIOHANSEL_FASTA } from '../../../../modules/nf-core/biohansel/main.nf'

workflow test_biohansel_paired_end {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    BIOHANSEL_PE ( input, [])
}

workflow test_biohansel_single_end {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
    ]

    BIOHANSEL_SE ( input, [])
}

workflow test_biohansel_assembly {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['genome']['genome_fasta_gz'], checkIfExists: true) ]
    ]

    BIOHANSEL_FASTA ( input, [])
}
