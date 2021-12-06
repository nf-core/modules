#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MASHTREE } from '../../../modules/mashtree/main.nf'

workflow test_mashtree {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
          file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
    ]

    MASHTREE ( input )
}
