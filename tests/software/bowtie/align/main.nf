#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE_BUILD } from '../../../../software/bowtie/build/main.nf' addParams( options: [:] )
include { BOWTIE_ALIGN } from '../../../../software/bowtie/align/main.nf' addParams( options: [:] )

workflow test_bowtie_align_single_end {

    def fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genomic.fasta", checkIfExists: true)
    BOWTIE_BUILD ( fasta )

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true) ] ]
    BOWTIE_ALIGN ( input, BOWTIE_BUILD.out.index )
}

workflow test_bowtie_align_paired_end {

    def fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genomic.fasta", checkIfExists: true)
    BOWTIE_BUILD ( fasta )

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_2.fastq.gz", checkIfExists: true) ] ]
    BOWTIE_ALIGN ( input, BOWTIE_BUILD.out.index )
}
