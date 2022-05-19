#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARIMA_2READSBAMCOMBINER } from '../../../../modules/arima/2readsbamcombiner/main.nf'

workflow test_arima_2readsbamcombiner {

    data = [
        [ id:'test'], // meta map
        [
           file(params.tol_test_data['test']['dImpGla2']['genomic_data']['hic_aligned_1_bam'] , checkIfExists: true),
           file(params.tol_test_data['test']['dImpGla2']['genomic_data']['hic_aligned_2_bam'] , checkIfExists: true),
        ],
    ]

    qscore=0

    ARIMA_2READSBAMCOMBINER ( data, qscore )
}
