#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_MERGE } from '../../../../../modules/nf-core/bbmap/bbmerge/main.nf'

workflow test_bbmap_merge_paired_fastq {

    input = [ [ id:'test', single_end:false ], // meta map
                [
                    file( params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                    file( params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
                ]
            ]

    BBMAP_MERGE ( input ) 
}


workflow test_bbmap_merge_single_fastq {

    input = [ [ id:'test', single_end:true ], // meta map
                [
                    file( params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
                ]
            ]

    BBMAP_MERGE ( input )
}