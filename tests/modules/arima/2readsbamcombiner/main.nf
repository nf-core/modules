#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARIMA_2READSBAMCOMBINER } from '../../../../modules/arima/2readsbamcombiner/main.nf'

workflow test_arima_2readsbamcombiner {

    data = [
        [ id:'test'], // meta map
        [
           file('https://tol-testdata.cog.sanger.ac.uk/test.single_end_1.bam', checkIfExists: true),
           file('https://tol-testdata.cog.sanger.ac.uk/test.single_end_2.bam', checkIfExists: true)
        ],
    ]

    qscore=0

    ARIMA_2READSBAMCOMBINER ( data, qscore )
}
