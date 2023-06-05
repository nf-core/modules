#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SHIGEIFINDER } from '../../../../modules/nf-core/shigeifinder/main.nf'
include { SHIGEIFINDER as SHIGEIFINDER_READS } from '../../../../modules/nf-core/shigeifinder/main.nf'

workflow test_shigeifinder_assembly {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    SHIGEIFINDER ( input )
}

workflow test_shigeifinder_reads {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
          file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
    ]

    SHIGEIFINDER_READS ( input )
}
