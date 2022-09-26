#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BIOBAMBAM_BAMSORMADUP } from '../../../../../modules/nf-core/biobambam/bamsormadup/main.nf'

workflow test_biobambam_bamsormadup_multi_input {

    input = [
        [ id:'test', single_end:false ], // meta map
        [file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true), file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)],
    ]

    BIOBAMBAM_BAMSORMADUP ( input, [] )
}

workflow test_biobambam_bamsormadup_single_input {

    input = [
        [ id:'test', single_end:false ], // meta map
        [file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)],
    ]

    BIOBAMBAM_BAMSORMADUP ( input, [] )
}
