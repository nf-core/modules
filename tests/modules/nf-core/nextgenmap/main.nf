#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NEXTGENMAP } from '../../../../modules/nf-core/nextgenmap/main.nf'

//
// Test with single-end data
//
workflow test_nextgenmap_single {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    NEXTGENMAP ( input, fasta )
}

//
// Test with paired-end data
//
workflow test_bwamem2_mem_paired_end {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    NEXTGENMAP ( input, fasta )
}
