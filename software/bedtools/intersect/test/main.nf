#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.intersect_args = '' //'-bed -c -f 0.20'

include check_output from  '../../../../tests/functions/check_process_outputs.nf' // params(params)
include INTERSECT_BED from '../main.nf' params(params)

// Define input channels
ch_input_1 = Channel.fromPath('./input_data/A.bed')
ch_input_2 = Channel.fromPath('./input_data/B.bed')

// Run the workflow
workflow {
    INTERSECT_BED(ch_input_1, ch_input_2, params.intersect_args)
    // .check_output()
}
