#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
include '../../../nf-core/module_testing/check_process_outputs.nf' params(params)
include '../main.nf' params(params)

// Define input channels
input = '../../../test-datasets/tools/bwa/index/input/reference.fasta'

// Run the workflow
workflow {
    samtools_faidx(input)
    // .check_output()
}
