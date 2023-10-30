#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STECFINDER } from '../../../../modules/nf-core/stecfinder/main.nf'
include { STECFINDER as STECFINDER_READS } from '../../../../modules/nf-core/stecfinder/main.nf'

workflow test_stecfinder_fasta {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    STECFINDER ( input )
}

workflow test_stecfinder_reads {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
          file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ],
    ]

    STECFINDER_READS ( input )
}
