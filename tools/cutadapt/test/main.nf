#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
include '../../../nf-core/module_testing/check_process_outputs.nf' params(params)
include '../main.nf' params(params)

// Define input channels
readPaths = [
  [ sample: 'SRR4238351', 
    R1: '../../../test-datasets/tools/cutadapt/input/SRR396636.sra_1.fastq', 
    R2: '../../../test-datasets/tools/cutadapt/input/SRR396636.sra_2.fastq'
  ]
]
Channel
  .from(readPaths)
  .map { row -> tuple( row.sample_name, file(row.R1.trim()), file(row.R2.trim()) ) }
  .set { ch_read_files }

// Run the workflow
workflow {
    cutadapt(ch_read_files)
    // .check_output()
}
