#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NGSBITS_SAMPLEGENDER } from '../../../../../modules/nf-core/ngsbits/samplegender/main.nf'

workflow test_ngsbits_samplegender {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    fasta = [
        [ id:'reference'], // meta map
        []
    ]

    fai = [
        [ id:'reference'], // meta map
        []
    ]

    NGSBITS_SAMPLEGENDER (
        input,
        fasta,
        fai,
        "sry"
    )
}
