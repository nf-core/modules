#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQKIT_REPLACE } from '../../../../modules/seqkit/replace/main.nf'
include { SEQKIT_REPLACE as SEQKIT_REPLACEUNCOMP } from '../../../../modules/seqkit/replace/main.nf'

workflow test_seqkit_replace {

    input = [ [ id:'test' ],   // meta map
              [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
            ]

    SEQKIT_REPLACE ( input )
}

workflow test_seqkit_replace_uncomp {

    input = [ [ id:'test' ],   // meta map
              [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
            ]

    SEQKIT_REPLACEUNCOMP ( input )
}
