#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RASUSA } from '../../../software/rasusa/main.nf' addParams( options: ['suffix':'_100X'])

workflow test_rasusa {
    def input = []
    input = [ [ id:'test', single_end:false], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
              ]
            ]

    depth_cutoff = 100
    genome_size = "1000000b"

    RASUSA ( input, depth_cutoff, genome_size )
}
