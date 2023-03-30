#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOURMASH_SKETCH } from '../../../../../modules/nf-core/sourmash/sketch/main.nf'
include { SOURMASH_INDEX } from '../../../../../modules/nf-core/sourmash/index/main.nf'

workflow test_sourmash_index {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    SOURMASH_SKETCH ( input )

    SOURMASH_INDEX ( 
        SOURMASH_SKETCH.out.signatures
    )
}
