#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQKIT_SPLIT2 as SEQKIT_SPLIT2_LENGTH } from '../../../../software/seqkit/split2/main.nf' addParams( options: ['args': '--by-length 26K'] )
include { SEQKIT_SPLIT2 as SEQKIT_SPLIT2_SIZE } from '../../../../software/seqkit/split2/main.nf' addParams( options: ['args': '--by-size 200' ] )
include { SEQKIT_SPLIT2 as SEQKIT_SPLIT2_PART } from '../../../../software/seqkit/split2/main.nf' addParams( options: ['args': '--by-part 3'] )

workflow test_seqkit_split2_single_end_length {
    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              file("${launchDir}/tests/data/fastq/dna/SRR396636_R1.fastq.gz", checkIfExists: true) ]

    SEQKIT_SPLIT2_LENGTH ( input )
}

workflow test_seqkit_split2_single_end_size {
    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              file("${launchDir}/tests/data/fastq/dna/SRR396636_R1.fastq.gz", checkIfExists: true) ]

    SEQKIT_SPLIT2_SIZE ( input )
}

workflow test_seqkit_split2_single_end_part {
    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              file("${launchDir}/tests/data/fastq/dna/SRR396636_R1.fastq.gz", checkIfExists: true) ]

    SEQKIT_SPLIT2_PART ( input )
}

workflow test_seqkit_split2_paired_end_length {
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/fastq/dna/SRR396636_*", checkIfExists: true) ]

    SEQKIT_SPLIT2_LENGTH ( input )
}

workflow test_seqkit_split2_paired_end_size {
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/fastq/dna/SRR396636_*", checkIfExists: true) ]

    SEQKIT_SPLIT2_SIZE ( input )
}

workflow test_seqkit_split2_paired_end_part {
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/fastq/dna/SRR396636_*", checkIfExists: true) ]

    SEQKIT_SPLIT2_PART ( input )
}