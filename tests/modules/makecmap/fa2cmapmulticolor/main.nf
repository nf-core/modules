#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAKECMAP_FA2CMAPMULTICOLOR } from '../../../../modules/makecmap/fa2cmapmulticolor/main.nf'

workflow test_makecmap_fa2cmapmulticolor {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    enzyme = "bspq1"

    MAKECMAP_FA2CMAPMULTICOLOR ( input, enzyme )
}
