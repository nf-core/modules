#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MASH_SKETCH } from '../../../../modules/mash/sketch/main.nf'

workflow test_mash_sketch {
    
    input = [ [ id:'test', single_end:false], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
              ]
            ]

    MASH_SKETCH ( input )
}
