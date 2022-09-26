#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KAT_HIST } from '../../../../../modules/nf-core/kat/hist/main.nf'

workflow test_kat_hist_single_end {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_1_fastq_gz'], checkIfExists: true)
    ]

    KAT_HIST ( input )
}

workflow test_kat_hist_paired_end {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test2_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_2_fastq_gz'], checkIfExists: true),
        ]
    ]

    KAT_HIST ( input )
}
