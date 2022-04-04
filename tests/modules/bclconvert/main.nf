#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCLCONVERT } from '../../../modules/bclconvert/main.nf'

workflow test_bclconvert {
    // TODO: add testdata or use stub
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    BCLCONVERT (samplesheet, run_dir )
}
