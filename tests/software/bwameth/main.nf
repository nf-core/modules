#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMETH_INDEX } from '../../../software/bwameth/index/main.nf' addParams( options: [:] )
include { BWAMETH_ALIGN as BWAMETH_ALIGN_SE } from '../../../software/bwameth/align/main.nf' addParams( options: [ publish_dir:'test_single_end' ] )
include { BWAMETH_ALIGN as BWAMETH_ALIGN_PE } from '../../../software/bwameth/align/main.nf' addParams( options: [ publish_dir:'test_paired_end' ] )

workflow test_bwameth_index {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true) ]

    BWAMETH_INDEX ( input )
}

/*
 * Test with single-end data
 */
workflow test_bwameth_align_single_end {

    def input = []
    def index = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/fastq/methylated_dna/Ecoli_10K_methylated_R1.fastq.gz", checkIfExists: true) ] ]
    index = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/index/E_coli/bwameth/NC_010473.fa.bwameth.c2t", checkIfExists: true),
                file("${launchDir}/tests/data/index/E_coli/bwameth/NC_010473.fa.bwameth.c2t.amb", checkIfExists: true),
                file("${launchDir}/tests/data/index/E_coli/bwameth/NC_010473.fa.bwameth.c2t.ann", checkIfExists: true),
                file("${launchDir}/tests/data/index/E_coli/bwameth/NC_010473.fa.bwameth.c2t.bwt", checkIfExists: true),
                file("${launchDir}/tests/data/index/E_coli/bwameth/NC_010473.fa.bwameth.c2t.sa", checkIfExists: true),
                file("${launchDir}/tests/data/index/E_coli/bwameth/NC_010473.fa.bwameth.c2t.pac", checkIfExists: true) ] ]

    BWAMETH_ALIGN_SE (
        input,
        index
    )
}

/*
 * Test with paired-end data
 */
workflow test_bwameth_align_paired_end {

    def input = []
    def index = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/fastq/methylated_dna/Ecoli_10K_methylated_R1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/fastq/methylated_dna/Ecoli_10K_methylated_R2.fastq.gz", checkIfExists: true) ] ]
    index = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/index/E_coli/bwameth/NC_010473.fa.bwameth.c2t", checkIfExists: true),
                file("${launchDir}/tests/data/index/E_coli/bwameth/NC_010473.fa.bwameth.c2t.amb", checkIfExists: true),
                file("${launchDir}/tests/data/index/E_coli/bwameth/NC_010473.fa.bwameth.c2t.ann", checkIfExists: true),
                file("${launchDir}/tests/data/index/E_coli/bwameth/NC_010473.fa.bwameth.c2t.bwt", checkIfExists: true),
                file("${launchDir}/tests/data/index/E_coli/bwameth/NC_010473.fa.bwameth.c2t.sa", checkIfExists: true),
                file("${launchDir}/tests/data/index/E_coli/bwameth/NC_010473.fa.bwameth.c2t.pac", checkIfExists: true) ] ]

    BWAMETH_ALIGN_PE (
        input,
        index
    )
}
