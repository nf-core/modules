#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.depth_cutoff = 100
params.genome_size = '1000000b'

include { RASUSA } from '../../../software/rasusa/main.nf' addParams( options: [:], depth_cutoff: params.depth_cutoff, genome_size:  params.genome_size)

workflow test_rasusa {
    def input = []
    input = [ [ id:'test', single_end:false], // meta map
              [file('https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/umi-dna/qiaseq/SRR7545951-small_1.fastq.gz', checkIfExists: true),
              file('https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/umi-dna/qiaseq/SRR7545951-small_2.fastq.gz', checkIfExists: true)
              ]
            ]

    RASUSA ( input )
}
