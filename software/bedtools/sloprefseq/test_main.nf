#!/usr/bin/env nextflow

nextflow.preview.dsl = 2
params.out_dir = "test_output"
params.upstream = 1
params.downstream = 10

include BEDTOOLS_SLOPEREFSEQ from '../main.nf' params(params)

// Define input channels
ch_input_1 = Channel.fromPath('./tests/data/bed/A.bed')
ch_input_2 = Channel.fromPath('./tests/data/bed/genome.sizes')

def additional_params_map = [:]

// Run the workflow
workflow {
    BEDTOOLS_SLOPEREFSEQ(ch_input_1, ch_input_2, additional_params_map)
}
