#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EMBOSS_REVSEQ } from '../../../../../modules/nf-core/emboss/revseq/main.nf'

workflow test_emboss_revseq {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    EMBOSS_REVSEQ ( input )
}
