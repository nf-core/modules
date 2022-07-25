#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTP } from '../../../modules/fastp/main.nf'

//
// Test with single-end data
//
workflow test_fastp_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    save_trimmed_fail = false
    save_merged       = false

    FASTP ( input, save_trimmed_fail, save_merged )
}

//
// Test with paired-end data
//
workflow test_fastp_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    save_trimmed_fail = false
    save_merged       = false

    FASTP ( input, save_trimmed_fail, save_merged )
}

//
// Test with single-end data with saving trimming fails
//
workflow test_fastp_single_end_trim_fail {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    save_trimmed_fail = true
    save_merged       = false

    FASTP ( input, save_trimmed_fail, save_merged )
}

//
// Test with paired-end data with saving trimming fails
//
workflow test_fastp_paired_end_trim_fail {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    save_trimmed_fail = true
    save_merged       = false

    FASTP ( input, save_trimmed_fail, save_merged )
}

//
// Test with paired-end data with merging
//
workflow test_fastp_paired_end_merged {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    save_trimmed_fail = false
    save_merged       = true

    FASTP ( input, save_trimmed_fail, save_merged )
}
