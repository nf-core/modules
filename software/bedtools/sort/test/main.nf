#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.out_dir = "test_output"
params.fastqc_args = ''
params.publish_dir_mode = "copy"
params.bedtools_sort_args = '' //'-sizeD'

include BEDTOOLS_SORT from '../main.nf' params(params)

// Define input channels
ch_input = Channel.fromPath('./input_data/A.bed')

// Run the workflow
workflow {
    BEDTOOLS_SORT(ch_input, params.bedtools_sort_args)
}
