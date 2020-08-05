#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

include check_output from  '../../../tests/functions/check_process_outputs.nf'
include tcoffee from '../main.nf'

// Define input channels
fasta = Channel.fromPath('../../../test-datasets/tools/tcoffee/input/BBA0001.tfa')

// Run the workflow
workflow {
    tcoffee(fasta)
    // .check_output()
}
