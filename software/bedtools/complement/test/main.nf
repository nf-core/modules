#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.bedtools_complement_args = ''

include check_output from  '../../../../tests/functions/check_process_outputs.nf' // params(params)
include BEDTOOLS_COMPLEMENT from '../main.nf' params(params)

// Define input channels
ch_input = Channel.fromPath('./input_data/A.bed')
chrom_sizes = Channel.fromPath('./input_data/genome.sizes')

// Run the workflow
workflow {
    BEDTOOLS_COMPLEMENT(ch_input, chrom_sizes, params.bedtools_complement_args)
    // .check_output()
}
