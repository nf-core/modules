#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//include { SEQKIT_SPLIT2_LENGTH } from '../../../../software/SEQKIT/SPLIT2/main.nf' addParams( options: ['args': '--by-length 26K'] )
//include { SEQKIT_SPLIT2_SIZE } from '../../../../software/SEQKIT/SPLIT2/main.nf' addParams( options: ['args': '--by-size 5000' ] )
include { SEQKIT_SPLIT2 as SEQKIT_SPLIT2_PART } from '../../../../software/SEQKIT/SPLIT2/main.nf' addParams( options: ['args': '--by-part 2'] )

// workflow test_seqkit_split2_length_single_end {
//     def input = []
//     input = [ [ id:'test', single_end:true ], // meta map
//               file("${launchDir}/tests/data/dna/SRR396636_R1.fastq.gz", checkIfExists: true) ]

//     SEQKIT_SPLIT2_LENGTH ( input )
// }

// workflow test_seqkit_split2_size_single_end {
//     def input = []
//     input = [ [ id:'test', single_end:true ], // meta map
//               file("${launchDir}/tests/data/dna/SRR396636_R1.fastq.gz", checkIfExists: true) ]

//     SEQKIT_SPLIT2_SIZE ( input )
// }

workflow test_seqkit_split2_part_single_end {
    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              file("${launchDir}/tests/data/dna/SRR396636_R1.fastq.gz", checkIfExists: true) ]

    SEQKIT_SPLIT2_PART ( input )
}

// workflow test_seqkit_split2_length_paired_end {
//     def input = []
//     input = [ [ id:'test', single_end:false ], // meta map
//               file("${launchDir}/tests/data/dna/SRR396636_*", checkIfExists: true) ]

//     SEQKIT_SPLIT2_LENGTH ( input )
// }

// workflow test_seqkit_split2_size_paired_end {
//     def input = []
//     input = [ [ id:'test', single_end:false ], // meta map
//               file("${launchDir}/tests/data/dna/SRR396636_*", checkIfExists: true) ]

//     SEQKIT_SPLIT2_SIZE ( input )
// }

// workflow test_seqkit_split2_part_paired_end {
//     def input = []
//     input = [ [ id:'test', single_end:false ], // meta map
//               file("${launchDir}/tests/data/dna/SRR396636_*", checkIfExists: true) ]

//     SEQKIT_SPLIT2_PART ( input )
// }