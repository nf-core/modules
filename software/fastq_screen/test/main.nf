#!/usr/bin/env nextflow
nextflow.preview.dsl = 2

params.outdir = "."
params.fastq_screen_args = ''
// fastq_screen_args are best passed in to the workflow in the following manner:
// --fastq_screen_args="--subset 200000 --force"

params.verbose = false

if (params.verbose){
    println ("[WORKFLOW] FASTQ SCREEN ARGS ARE: " + params.fastq_screen_args)
}

// TODO: include '../../../tests/functions/check_process_outputs.nf'
include '../main.nf'

// Define input channels

ch_read_files = Channel 
  .fromFilePairs('../../../test-datasets/Ecoli*{1,2}.fastq.gz',size:-1)
  // .view()  // to check whether the input channel works

// Run the workflow
workflow {
    main:
        FASTQ_SCREEN(ch_read_files, params.outdir, params.fastq_screen_args, params.verbose)

    // TODO .check_output()
}
