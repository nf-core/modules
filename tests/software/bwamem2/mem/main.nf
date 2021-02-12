#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMEM2_MEM } from '../../../../software/bwamem2/mem/main.nf' addParams( options: [:] )

/*
 * Test with single-end data
 */
workflow test_bwamem2_mem_single_end {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/fastq/dna/Ecoli_DNA_R1.fastq.gz", checkIfExists: true) ] ]

    BWAMEM2_MEM (
        input,
        file("${launchDir}/tests/data/index/E_coli/bwamem2/", checkIfExists: true)
    )
}

/*
 * Test with paired-end data
 */
workflow test_bwamem2_mem_paired_end {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/fastq/dna/Ecoli_DNA_R1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/fastq/dna/Ecoli_DNA_R2.fastq.gz", checkIfExists: true) ] ]

    BWAMEM2_MEM (
        input,
        file("${launchDir}/tests/data/index/E_coli/bwamem2/", checkIfExists: true)
    )
}
