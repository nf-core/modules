#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.complementbed_args = ''

include check_output from  '../../../../tests/functions/check_process_outputs.nf' // params(params)
include COMPLEMENT_BED from '../main.nf' params(params)

// Define input channels
ch_input = Channel.fromPath('./input_data/A.bed')
chrom_sizes = Channel.fromPath('./input_data/genome.sizes')

// Run the workflow
workflow {
    COMPLEMENT_BED(ch_input, chrom_sizes, params.complementbed_args)
    // .check_output()
}
