#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RASUSA } from '../../../modules/rasusa/main.nf'

workflow test_rasusa {
    input = [ [ id:'test', single_end:false], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
              ],
              "1000000b"
            ]

    depth_cutoff = 100

    RASUSA ( input, depth_cutoff )
}
