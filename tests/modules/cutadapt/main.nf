#!/usr/bin/env nextflow



include { CUTADAPT } from '../../../modules/cutadapt/main.nf'

//
// Test with single-end data
//
workflow test_cutadapt_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]

    CUTADAPT ( input )
}

//
// Test with paired-end data
//
workflow test_cutadapt_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    CUTADAPT ( input )
}

