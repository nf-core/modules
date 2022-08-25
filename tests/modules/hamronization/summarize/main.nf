#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HAMRONIZATION_DEEPARG                                 } from '../../../../modules/hamronization/deeparg/main.nf'
include { HAMRONIZATION_DEEPARG as HAMRONIZATION_DEEPARG_SECOND } from '../../../../modules/hamronization/deeparg/main.nf'
include { HAMRONIZATION_SUMMARIZE                               } from '../../../../modules/hamronization/summarize/main.nf'

workflow test_hamronization_summarize {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_mapping_potential_arg'], checkIfExists: true),
    ]

    input2 = [
        [ id:'test2', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_mapping_potential_arg'], checkIfExists: true),
    ]

    HAMRONIZATION_DEEPARG ( input, 'tsv', '1.0.2', '2'  )
    HAMRONIZATION_DEEPARG_SECOND ( input2, 'tsv', '1.0.2', '2' )

    ch_deeparg_run_one = HAMRONIZATION_DEEPARG.out.tsv
    ch_deeparg_run_two = HAMRONIZATION_DEEPARG_SECOND.out.tsv

    ch_deeparg_run_one
        .mix( ch_deeparg_run_two )
        .map{
            [ it[1] ]
        }
        .collect()
        .set { ch_input_for_summarize }

    HAMRONIZATION_SUMMARIZE ( ch_input_for_summarize , 'json' )
}
