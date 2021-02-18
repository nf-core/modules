#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMETH_ALIGN as BWAMETH_ALIGN_SE } from '../../../../software/bwameth/align/main.nf' addParams( options: [ publish_dir:'test_single_end' ] )
include { BWAMETH_ALIGN as BWAMETH_ALIGN_PE } from '../../../../software/bwameth/align/main.nf' addParams( options: [ publish_dir:'test_paired_end' ] )

/*
 * Test with single-end data
 */
workflow test_bwameth_align_single_end {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/fastq/methylated_dna/Ecoli_10K_methylated_R1.fastq.gz", checkIfExists: true) ] ]

    BWAMETH_ALIGN_SE (
        input,
        file("${launchDir}/tests/data/index/E_coli/bwameth", checkIfExists: true)
    )
}

/*
 * Test with paired-end data
 */
workflow test_bwameth_align_paired_end {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/fastq/methylated_dna/Ecoli_10K_methylated_R1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/fastq/methylated_dna/Ecoli_10K_methylated_R2.fastq.gz", checkIfExists: true) ] ]

    BWAMETH_ALIGN_PE (
        input,
        file("${launchDir}/tests/data/index/E_coli/bwameth", checkIfExists: true)
    )
}
