#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STADENIOLIB_CRAMFILTER } from '../../../../../modules/nf-core/stadeniolib/cramfilter/main.nf'

workflow test_stadeniolib_cramfilter {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
    ]

    STADENIOLIB_CRAMFILTER ( input, 1, 1000 )
}