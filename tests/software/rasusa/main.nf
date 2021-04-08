#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RASUSA } from '../../../software/rasusa/main.nf' addParams( options: [:])

workflow test_rasusa {
    def input = []
    input = [ [ id:'test', single_end:false], // meta map
              [file('https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/umi-dna/qiaseq/SRR7545951-small_1.fastq.gz', checkIfExists: true),
              file('https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/umi-dna/qiaseq/SRR7545951-small_2.fastq.gz', checkIfExists: true)
              ]
            ]

    depth_cutoff = 100

    genome_size = "1000000b"

    RASUSA ( input, depth_cutoff, genome_size )
}
