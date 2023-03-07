#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SENTIEON_BWAINDEX } from '../../../../../modules/nf-core/sentieon/bwaindex/main.nf'
include { SENTIEON_BWAMEM } from '../../../../../modules/nf-core/sentieon/bwamem/main.nf'

//
// Test with single-end data
//
workflow test_sentieon_bwamem_single_end {
    input_ch = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta_file = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai_file = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    fasta_ch = [
        [ id: 'test' ],
        fasta_file
    ]

    SENTIEON_BWAINDEX ( fasta_ch )
    SENTIEON_BWAMEM ( input_ch, SENTIEON_BWAINDEX.out.index, fasta_file, fasta_fai_file )
}

/*
//
// Test with single-end data and sort
//
workflow test_bwa_mem_single_end_sort {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BWA_INDEX ( fasta )
    BWA_MEM ( input, BWA_INDEX.out.index, true )
}

//
// Test with paired-end data
//
workflow test_bwa_mem_paired_end {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BWA_INDEX ( fasta )
    BWA_MEM ( input, BWA_INDEX.out.index, false )
}

//
// Test with paired-end data and sort
//
workflow test_bwa_mem_paired_end_sort {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BWA_INDEX ( fasta )
    BWA_MEM ( input, BWA_INDEX.out.index, true )
}
*/
