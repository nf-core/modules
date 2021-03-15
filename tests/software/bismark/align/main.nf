#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOME_PREPARATION } from '../../../../software/bismark/genome_preparation/main.nf' addParams( options: [:] )
include { BISMARK_ALIGN as BISMARK_ALIGN_SE } from '../../../../software/bismark/align/main.nf'  addParams( options: [ publish_dir:'test_single_end' ] )
include { BISMARK_ALIGN as BISMARK_ALIGN_PE } from '../../../../software/bismark/align/main.nf'  addParams( options: [ publish_dir:'test_paired_end' ] )

/*
 * Test with single-end data
 */
workflow test_bismark_align_single_end {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_methylated_1.fastq.gz", checkIfExists: true) ] ]


    BISMARK_GENOME_PREPARATION ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) )
    BISMARK_ALIGN_SE (
        input,
        BISMARK_GENOME_PREPARATION.out.index
    )
}

/*
 * Test with paired-end data
 */
workflow test_bismark_align_paired_end {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_methylated_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_methylated_2.fastq.gz", checkIfExists: true) ] ]


    BISMARK_GENOME_PREPARATION ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) )
    BISMARK_ALIGN_PE (
        input,
        BISMARK_GENOME_PREPARATION.out.index
    )
}
