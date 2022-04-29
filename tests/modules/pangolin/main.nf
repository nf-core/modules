#!/usr/bin/env nextflow



include { PANGOLIN } from '../../../modules/pangolin/main.nf'

workflow test_pangolin {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ] ]

    PANGOLIN ( input )
}
