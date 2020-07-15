#!/usr/bin/env nextflow
nextflow.preview.dsl = 2

params.out_dir = "test_output"
params.fastqc_args = ''
params.publish_dir_mode = "copy"

include { FASTQC } from '../main.nf'

/**
 * Test if FASTQC runs with single-end data
 */
workflow test_single_end {
    input_files = Channel.fromPath("data/test_single_end.fastq.gz")
                    .map {f -> [f.baseName, f]}
    FASTQC(input_files)
}

/**
 * Test if FASTQC runs with paired end data
 */
workflow test_paired_end {
    input_files = Channel.fromFilePairs("data/test_R{1,2}.fastq.gz")
    FASTQC(input_files)
}

workflow {
    test_single_end()
    test_paired_end()
}
