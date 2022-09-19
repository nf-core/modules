#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ABRICATE_RUN } from '../../../../modules/abricate/run/main.nf'
include { ABRICATE_SUMMARY } from '../../../../modules/abricate/summary/main.nf'

workflow test_abricate_summary {

    inputs = [
        tuple([ id:'test1', single_end:false ], // meta map
              file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)),
        tuple([ id:'test2', single_end:false ],
              file(params.test_data['haemophilus_influenzae']['genome']['genome_fna_gz'], checkIfExists: true))
    ]

    ABRICATE_RUN ( Channel.fromList(inputs) )
    ABRICATE_SUMMARY ( 
        ABRICATE_RUN.out.report.collect{ meta, report -> report }.map{ report -> [[ id: 'test_summary'], report]}
    )
}
