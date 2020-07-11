#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.outdir = "."             // gets set in the nextflow.config files (to './results/trim_galore')
params.verbose = false
params.trim_galore_args = ''
// trim_galore_args are best passed into the workflow in the following manner, e.g.:
// --trim_galore_args="--clip_r1 10 --clip_r2 15 -j 2"

if (params.verbose){
    println ("[WORKFLOW] TRIM GALORE ARGS: "      + params.trim_galore_args)
}

// TODO: check the output files in some way
// include '../../../tests/functions/check_process_outputs.nf'
include '../main.nf'  // params (clip_r1: 6, clip_r2: 10) // how to pass additional parameters

ch_read_files = Channel
  .fromFilePairs('../../../test-datasets/test*{1,2}.fastq.gz',size:-1)
  // .view()  // to check whether the input channel works

workflow {

    main:
        TRIM_GALORE (ch_read_files, params.outdir, params.trim_galore_args, params.verbose)

}
