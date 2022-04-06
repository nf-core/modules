#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCLCONVERT } from '../../../modules/bclconvert/main.nf'

process STUB_BCLCONVERT_INPUT {
    output:
    path "SampleSheet.csv"          ,emit: samplesheet
    path "DDMMYY_SERIAL_FLOWCELL"   ,emit: run_dir

    stub:
    """
    mkdir DDMMYY_SERIAL_FLOWCELL
    echo "SampleSheet" > SampleSheet.csv
    """
}

workflow test_bclconvert {
    STUB_BCLCONVERT_INPUT ()
    BCLCONVERT (STUB_BCLCONVERT_INPUT.out.samplesheet, STUB_BCLCONVERT_INPUT.out.run_dir)
}
