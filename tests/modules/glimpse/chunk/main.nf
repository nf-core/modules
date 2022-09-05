#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GLIMPSE_CHUNK } from '../../../../modules/glimpse/chunk/main.nf'

workflow test_glimpse_chunk {
    input = [
        [ id:'test', single_end:false ], // meta map
        [file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz_tbi'], checkIfExists: true)]
    ]
    GLIMPSE_CHUNK (input,"chr21")
}
