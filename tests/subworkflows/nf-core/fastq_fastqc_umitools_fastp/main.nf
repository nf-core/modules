#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQ_FASTQC_UMITOOLS_FASTP } from '../../../../subworkflows/nf-core/fastq_fastqc_umitools_fastp/main.nf'

workflow test_fastq_fastqc_umitools_trimgalore_paired {
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
    umi_discard_read  = 1
    skip_trimming     = false
    adapter_fasta     = []
    save_trimmed_fail = false
    save_merged       = false
    min_trimmed_reads = 1

    FASTQ_FASTQC_UMITOOLS_FASTP (
        input,
        skip_fastqc,
        with_umi,
        skip_umi_extract,
        umi_discard_read,
        skip_trimming,
        adapter_fasta,
        save_trimmed_fail,
        save_merged,
        min_trimmed_reads
    )
}
