#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { RGI_MAIN } from '../../../../../modules/nf-core/rgi/main/main.nf'
include { HAMRONIZATION_RGI } from '../../../../../modules/nf-core/hamronization/rgi/main.nf'

workflow test_hamronization_rgi {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['haemophilus_influenzae']['genome']['genome_fna_gz'], checkIfExists: true)
    ]

    RGI_MAIN ( input )
    HAMRONIZATION_RGI ( RGI_MAIN.out.tsv, 'tsv', '1.0.2', '3.2.3' )
}
