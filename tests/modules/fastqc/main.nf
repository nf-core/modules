#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC } from '../../../modules/fastqc/main.nf'

//
// Test with single-end data
//
workflow test_fastqc_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]

    FASTQC ( input )
}

//
// Test with paired-end data
//
workflow test_fastqc_paired_end {
    input = [ [id: 'test', single_end: false], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    FASTQC ( input )
}
