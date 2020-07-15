#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.out_dir = "test_output"
params.fastqc_args = ''
params.publish_dir_mode = "copy"
params.bedtools_complement_args = ''

include BEDTOOLS_COMPLEMENT from '../main.nf' params(params)

// Define input channels
ch_input = Channel.fromPath('./input_data/A.bed')
chrom_sizes = Channel.fromPath('./input_data/genome.sizes')

// Run the workflow
workflow {
    BEDTOOLS_COMPLEMENT(ch_input, chrom_sizes, params.bedtools_complement_args)
}
