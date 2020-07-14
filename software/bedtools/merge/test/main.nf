#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.bedtools_merge_args = '' //''-s -c 6 -o distinct'

include check_output from  '../../../../tests/functions/check_process_outputs.nf' // params(params)
include BEDTOOLS_MERGE from '../main.nf' params(params)

// Define input channels
ch_input = Channel.fromPath('../../test_data_set/A.bed')
//ch_input = Channel.fromPath('./input_data/JK2067_downsampled_s0.1.bam')

// Run the workflow
workflow {
    BEDTOOLS_MERGE(ch_input, params.bedtools_merge_args)
    // .check_output()
}
