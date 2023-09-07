#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALIGN_ALIGN } from '../../../../../modules/nf-core/kalign/align/main.nf'

workflow test_kalign_align {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]

    KALIGN_ALIGN ( input )
}