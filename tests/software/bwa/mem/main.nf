#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX } from '../../../../software/bwa/index/main.nf' addParams( options: [:] )
include { BWA_MEM } from '../../../../software/bwa/mem/main.nf' addParams( options: [:] )

/*
 * Test with single-end data
 */
workflow test_bwa_mem_single_end {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true) ] ]

    BWA_INDEX ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) )
    BWA_MEM (
        input,
        BWA_INDEX.out.index
    )
}

/*
 * Test with paired-end data
 */
workflow test_bwa_mem_paired_end {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_2.fastq.gz", checkIfExists: true) ] ]

    BWA_INDEX ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) )
    BWA_MEM (
        input,
        BWA_INDEX.out.index
    )
}
