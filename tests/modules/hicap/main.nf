#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HICAP } from '../../../modules/hicap/main.nf'

workflow test_hicap {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['haemophilus_influenzae']['genome']['genome_fna_gz'], checkIfExists: true)
    ]
    database_dir = []
    model_fp = []

    HICAP ( input, database_dir, model_fp )
}
