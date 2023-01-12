#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_SENDSKETCH } from '../../../../../modules/nf-core/bbmap/sendsketch/main.nf'

workflow test_bbmap_sendsketch {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['genome_fasta'], checkIfExists: true)
    ]

    BBMAP_SENDSKETCH ( input )
}
