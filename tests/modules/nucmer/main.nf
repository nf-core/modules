#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NUCMER } from '../../../modules/nucmer/main.nf'

workflow test_nucmer {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true) ]

    NUCMER ( input )
}
