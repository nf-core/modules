#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.bedtools_sort_args = '' //'-sizeD'

include check_output from  '../../../../tests/functions/check_process_outputs.nf' // params(params)
include BEDTOOLS_SORT from '../main.nf' params(params)

// Define input channels
ch_input = Channel.fromPath('./input_data/A.bed')

// Run the workflow
workflow {
    BEDTOOLS_SORT(ch_input, params.bedtools_sort_args)
    // .check_output()
}
