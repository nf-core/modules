#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMEM2_INDEX } from '../../../../modules/bwamem2/index/main.nf'
include { BWAMEM2_MEM   } from '../../../../modules/bwamem2/mem/main.nf'

//
// Test with single-end data
//
workflow test_bwamem2_mem_single_end {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BWAMEM2_INDEX ( fasta )
    BWAMEM2_MEM ( input, BWAMEM2_INDEX.out.index, false )
}

//
// Test with single-end data and sort
//
workflow test_bwamem2_mem_single_end_sort {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BWAMEM2_INDEX ( fasta )
    BWAMEM2_MEM ( input, BWAMEM2_INDEX.out.index, true )
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

    BWAMEM2_INDEX ( fasta )
    BWAMEM2_MEM ( input, BWAMEM2_INDEX.out.index, false )
}

//
// Test with paired-end data and sort
//
workflow test_bwamem2_mem_paired_end_sort {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BWAMEM2_INDEX ( fasta )
    BWAMEM2_MEM ( input, BWAMEM2_INDEX.out.index, true )
}
