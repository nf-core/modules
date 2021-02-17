#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQKIT_SPLIT2 } from '../../../../software/SEQKIT/SPLIT2/main.nf' addParams( options: [:] )

workflow test_seqkit_split2_single_end {
    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              file("${launchDir}/tests/data/dna/SRR396636_R1.fastq.gz", checkIfExists: true) ]

    SEQKIT_SPLIT2 ( input )
}

workflow test_seqkit_split2_paired_end {
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/dna/SRR396636_*", checkIfExists: true) ]

    SEQKIT_SPLIT2 ( input )
}