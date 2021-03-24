#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMETH_INDEX                     } from '../../../../software/bwameth/index/main.nf' addParams( options: [:]                               )
include { BWAMETH_ALIGN as BWAMETH_ALIGN_SE } from '../../../../software/bwameth/align/main.nf' addParams( options: [ publish_dir:'test_single_end' ] )
include { BWAMETH_ALIGN as BWAMETH_ALIGN_PE } from '../../../../software/bwameth/align/main.nf' addParams( options: [ publish_dir:'test_paired_end' ] )

/*
 * Test with single-end data
 */
workflow test_bwameth_align_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_methylated_1.fastq.gz", checkIfExists: true) ] 
            ]
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)

    BWAMETH_INDEX ( fasta )
    BWAMETH_ALIGN_SE ( input, BWAMETH_INDEX.out.index )
}

/*
 * Test with paired-end data
 */
workflow test_bwameth_align_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_methylated_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/fastq/test_methylated_2.fastq.gz", checkIfExists: true) ] 
            ]
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)

    BWAMETH_INDEX ( fasta )
    BWAMETH_ALIGN_PE ( input, BWAMETH_INDEX.out.index )
}
