#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from '../../../../software/bowtie2/build/main.nf' addParams( options: [:] )
include { BOWTIE2_ALIGN } from '../../../../software/bowtie2/align/main.nf' addParams( options: [:] )

workflow test_bowtie2_align_single_end {

    def fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)
    BOWTIE2_BUILD ( fasta )

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true) ] ]
    BOWTIE2_ALIGN ( input, BOWTIE2_BUILD.out.index )
}

workflow test_bowtie2_align_paired_end {

    def fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)
    BOWTIE2_BUILD ( fasta )

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_2.fastq.gz", checkIfExists: true) ] ]
    BOWTIE2_ALIGN ( input, BOWTIE2_BUILD.out.index )
}
