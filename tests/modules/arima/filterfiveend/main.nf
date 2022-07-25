#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARIMA_FILTERFIVEEND } from '../../../../modules/arima/filterfiveend/main.nf'

workflow test_arima_filterfiveend_single_end {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)
    ]

    ARIMA_FILTERFIVEEND ( input )
}

workflow test_arima_filterfiveend_paired_end {

    input = [
        [ id:'test', single_end:false ], // meta map
        [file(params.test_data['bacteroides_fragilis']['illumina']['test1_paired_end_bam'], checkIfExists: true),
        file(params.test_data['bacteroides_fragilis']['illumina']['test2_paired_end_bam'], checkIfExists: true)]
    ]

    ARIMA_FILTERFIVEEND ( input )
}
