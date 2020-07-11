#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
include '../main.nf' params(params)

// Define input channels

Channel
  .fromFilePairs('../../../test-datasets/tools/cutadapt/input/*_{1,2}.fastq' )
  .set { ch_read_files }

// Run the workflow
workflow {
    cutadapt(ch_read_files)
}
