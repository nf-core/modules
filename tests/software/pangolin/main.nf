#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PANGOLIN } from '../../../software/pangolin/main.nf' addParams( options: [:] )

workflow test_pangolin {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ] ]

    PANGOLIN ( input )
}
