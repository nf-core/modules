#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINIMAP2_ALIGN } from '../../../../software/minimap2/align/main.nf' addParams( options: [:] )

workflow test_minimap2_align_single_end {

    def fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true) ] ]
    MINIMAP2_ALIGN ( input, fasta )
}

workflow test_minimap2_align_paired_end {

    def fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_2.fastq.gz", checkIfExists: true) ] ]
    MINIMAP2_ALIGN ( input, fasta )
}
