#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPARG_DOWNLOADDATA } from '../../../../modules/deeparg/downloaddata/main.nf'
include { DEEPARG_PREDICT      } from '../../../../modules/deeparg/predict/main.nf'

workflow test_deeparg_predict {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true),
        'LS'
    ]

    DEEPARG_DOWNLOADDATA( )
    DEEPARG_PREDICT ( input, DEEPARG_DOWNLOADDATA.out.db )

}
