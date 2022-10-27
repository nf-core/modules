#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { FASTQC } from "$moduleDir/modules/nf-core/fastqc/main.nf"

//
// Test with single-end data
//
workflow test_fastqc_single_end {
    input = [
                [ id:'test', single_end:true ], // meta map
                [
                    file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
                ]
            ]

    FASTQC ( input )
}

//
// Test with paired-end data
//
workflow test_fastqc_paired_end {
    input = [
                [id: 'test', single_end: false], // meta map
                [
                    file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                    file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
                ]
            ]

    FASTQC ( input )
}

//
// Test with interleaved data
//
workflow test_fastqc_interleaved {
    input = [
                [id: 'test', single_end: false], // meta map
                [
                    file(params.test_data['sarscov2']['illumina']['test_interleaved_fastq_gz'], checkIfExists: true),
                ]
            ]

    FASTQC ( input )
}

//
// Test with bam data
//
workflow test_fastqc_bam {
    input = [
                [id: 'test', single_end: false], // meta map
                [
                    file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                ]
            ]

    FASTQC ( input )
}

//
// Test with multiple samples
//
workflow test_fastqc_multiple {
    input = [
                [id: 'test', single_end: false], // meta map
                [
                    file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                    file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true),
                    file(params.test_data['sarscov2']['illumina']['test2_1_fastq_gz'], checkIfExists: true),
                    file(params.test_data['sarscov2']['illumina']['test2_2_fastq_gz'], checkIfExists: true),

                ]
            ]

    FASTQC ( input )
}
