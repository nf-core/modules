#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../../../../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main.nf'

//
// Test with single-end data
//
workflow test_fastq_fastqc_umitools_trimgalore_single {
    input = [
                [ id:'test', single_end:true ], // meta map
                [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    skip_fastqc       = false
    with_umi          = true
    skip_umi_extract  = false
    skip_trimming     = false
    umi_discard_read  = 1
    min_trimmed_reads = 1

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE ( input, skip_fastqc, with_umi, skip_umi_extract, skip_trimming, umi_discard_read, min_trimmed_reads)
}

workflow test_fastq_fastqc_umitools_trimgalore_paired {
    input = [
                [ id:'test', single_end:false ], // meta map
                [
                    file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                    file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
                ]
            ]
    skip_fastqc       = false
    with_umi          = true
    skip_umi_extract  = false
    skip_trimming     = false
    umi_discard_read  = 1
    min_trimmed_reads = 1

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE ( input, skip_fastqc, with_umi, skip_umi_extract, skip_trimming, umi_discard_read, min_trimmed_reads)
}

workflow test_fastq_fastqc_umitools_trimgalore_ {
    input = [
                [ id:'test', single_end:false ], // meta map
                [
                    file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                    file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
                ]
            ]
    skip_fastqc       = false
    with_umi          = false
    skip_umi_extract  = false
    skip_trimming     = false
    umi_discard_read  = 1
    min_trimmed_reads = 1

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE ( input, skip_fastqc, with_umi, skip_umi_extract, skip_trimming, umi_discard_read, min_trimmed_reads)
}

workflow test_fastq_fastqc_umitools_trimgalore_umi {
    input = [
                [ id:'test', single_end:false ], // meta map
                [
                    file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                    file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
                ]
            ]
    skip_fastqc       = false
    with_umi          = false
    skip_umi_extract  = false
    skip_trimming     = false
    umi_discard_read  = 1
    min_trimmed_reads = 1

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE ( input, skip_fastqc, with_umi, skip_umi_extract, skip_trimming, umi_discard_read, min_trimmed_reads)
}
workflow test_fastq_fastqc_umitools_trimgalore_skip {
    input = [
                [ id:'test', single_end:false ], // meta map
                [
                    file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                    file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
                ]
            ]
    skip_fastqc       = true
    with_umi          = true
    skip_umi_extract  = true
    skip_trimming     = true
    umi_discard_read  = 0
    min_trimmed_reads = 1

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE ( input, skip_fastqc, with_umi, skip_umi_extract, skip_trimming, umi_discard_read, min_trimmed_reads)
}
