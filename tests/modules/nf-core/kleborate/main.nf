#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KLEBORATE } from '../../../modules/kleborate/main.nf'

workflow test_kleborate {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
          file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true),
          file(params.test_data['sarscov2']['illumina']['scaffolds_fasta'], checkIfExists: true)
        ]
    ]

    KLEBORATE ( input )
}
