#!/usr/bin/env nextflow
nextflow.preview.dsl = 2

params.outdir = "."             // gets set in nextflow.config file (as './results/fastqc')
params.fastqc_args = ''
params.verbose = false

// TODO: check the output files in some way
// include '../../../nf-core/module_testing/check_process_outputs.nf'
include '../main.nf'

// Define input channels
ch_read_files = Channel 
    .fromFilePairs('../../../test-datasets/test*{1,2}.fastq.gz',size:-1)
    // .view()  // to check whether the input channel works

// Run the workflow
workflow {
    FASTQC (ch_read_files, params.outdir, params.fastqc_args, params.verbose)
    // .check_output()
}
