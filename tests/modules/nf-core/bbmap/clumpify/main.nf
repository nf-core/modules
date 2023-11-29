#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_CLUMPIFY } from '../../../../../modules/nf-core/bbmap/clumpify/main.nf'

workflow test_bbmap_clumpify_single_end {

    input = [ [ id:'test', single_end:true ], // meta map
              [  file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]

    BBMAP_CLUMPIFY ( input )
}

workflow test_bbmap_clumpify_paired_end {

    input = [ [ id:'test', single_end:false ], // meta map
              [  file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                 file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    BBMAP_CLUMPIFY ( input )
}
