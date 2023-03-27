#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQ_TRIM_FASTP_FASTQC } from '../../../../subworkflows/nf-core/fastq_trim_fastp_fastqc/main.nf'

//
// Test with single-end data
//
workflow test_fastq_trim_fastp_fastqc_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    save_trimmed_fail = false
    save_merged       = false
    skip_fastqc       = false
    skip_fastp        = false

    FASTQ_TRIM_FASTP_FASTQC ( input, [], save_trimmed_fail, save_merged, skip_fastp, skip_fastqc )
}

//
// Test with paired-end data
//
workflow test_fastq_trim_fastp_fastqc_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    save_trimmed_fail = false
    save_merged       = false
    skip_fastqc       = false
    skip_fastp        = false

    FASTQ_TRIM_FASTP_FASTQC ( input, [], save_trimmed_fail, save_merged, skip_fastp, skip_fastqc )
}

//
// Test with intereleaved data
//
workflow test_fastq_trim_fastp_fastqc_interleaved {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_interleaved_fastq_gz'], checkIfExists: true) ]
            ]
    save_trimmed_fail = false
    save_merged       = false
    skip_fastqc       = false
    skip_fastp        = false

    FASTQ_TRIM_FASTP_FASTQC ( input, [], save_trimmed_fail, save_merged, skip_fastp, skip_fastqc )
}

//
// Test with single-end data with saving trimming fails
//
workflow test_fastq_trim_fastp_fastqc_single_end_trim_fail {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    save_trimmed_fail = true
    save_merged       = false
    skip_fastqc       = false
    skip_fastp        = false

    FASTQ_TRIM_FASTP_FASTQC ( input, [], save_trimmed_fail, save_merged, skip_fastp, skip_fastqc )
}

//
// Test with paired-end data with saving trimming fails
//
workflow test_fastq_trim_fastp_fastqc_paired_end_trim_fail {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    save_trimmed_fail = true
    save_merged       = false
    skip_fastqc       = false
    skip_fastp        = false

    FASTQ_TRIM_FASTP_FASTQC ( input, [], save_trimmed_fail, save_merged, skip_fastp, skip_fastqc )
}

//
// Test with paired-end data with merging
//
workflow test_fastq_trim_fastp_fastqc_paired_end_merged {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    save_trimmed_fail = false
    save_merged       = true
    skip_fastqc       = false
    skip_fastp        = false

    FASTQ_TRIM_FASTP_FASTQC ( input, [], save_trimmed_fail, save_merged, skip_fastp, skip_fastqc )
}

//
// Test with paired-end data with predefined adapter list
//
workflow test_fastq_trim_fastp_fastqc_paired_end_merged_adapterlist {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    adapter_fasta = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/fastp/adapters.fasta", checkIfExists: true)
    save_trimmed_fail = false
    save_merged       = true
    skip_fastqc       = false
    skip_fastp        = false

    FASTQ_TRIM_FASTP_FASTQC ( input, adapter_fasta, save_trimmed_fail, save_merged, skip_fastp, skip_fastqc )
}

//
// Test with paired-end data with skipping fastqc
//
workflow test_fastq_trim_fastp_fastqc_paired_end_skip_fastqc {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    save_trimmed_fail = false
    save_merged       = false
    skip_fastqc       = true
    skip_fastp        = false

    FASTQ_TRIM_FASTP_FASTQC ( input, [], save_trimmed_fail, save_merged, skip_fastp, skip_fastqc )
}

//
// Test with paired-end data with skipping fastp
//
workflow test_fastq_trim_fastp_fastqc_paired_end_skip_fastp {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    save_trimmed_fail = false
    save_merged       = false
    skip_fastqc       = false
    skip_fastp        = true

    FASTQ_TRIM_FASTP_FASTQC ( input, [], save_trimmed_fail, save_merged, skip_fastp, skip_fastqc )
}
