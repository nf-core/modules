#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PRODIGAL } from '../../../modules/prodigal/main.nf'

workflow test_prodigal {
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) 
            ]

    PRODIGAL ( input , "gff")
}
