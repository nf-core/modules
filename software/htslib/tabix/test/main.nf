#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
include '../../../tests/functions/check_process_outputs.nf' params(params)
include '../main.nf' params(params)

// Define input channels
input = '../../../test-datasets/tools/file.vcf.gz'

// Run the workflow
workflow {
    tabix_index(ch_read_files)
    // .check_output()
}
