#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.out_dir = "test_output"
params.fastqc_args = ''
params.publish_dir_mode = "copy"
params.bedtools_merge_args = '' //''-s -c 6 -o distinct'

include BEDTOOLS_MERGE from '../main.nf' params(params)

// Define input channels
ch_input = Channel.fromPath('./input_data/A.bed')
//ch_input = Channel.fromPath('./input_data/JK2067_downsampled_s0.1.bam')

// Run the workflow
workflow {
    BEDTOOLS_MERGE(ch_input, params.bedtools_merge_args)
}
