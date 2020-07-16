#!/usr/bin/env nextflow

nextflow.preview.dsl = 2
params.out_dir = "test_output"
params.fastqc_args = ''
params.publish_dir_mode = "copy"
params.intersect_args = '' //'-bed -c -f 0.20'

include check_output from  '../../../../tests/functions/check_process_outputs.nf' // params(params)
include INTERSECT_BED from '../main.nf' params(params)

// Define input channels
ch_input_1 = Channel.fromPath('./input_data/A.bed')
ch_input_2 = Channel.fromPath('./input_data/B.bed')

def additional_params_map = [:]

additional_params_map =  [ s: "",
                           f: 0.9 ]

// Run the workflow
workflow {
    INTERSECT_BED(ch_input_1, ch_input_2, additional_params_map)
}
