#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STRELKA_GERMLINE } from '../../../../software/strelka/germline/main.nf' addParams( options: [:] )

workflow test_strelka_germline {
    
    def input = []

    def fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)
    def fai   = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta.fai", checkIfExists: true)
    def target_bed = []
    input = [ [ id:'test'], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam.bai", checkIfExists: true) ]

    STRELKA_GERMLINE ( input, fasta, fai, target_bed )
}
workflow test_strelka_germline_target_bed {
    
    def input = []

    def fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)
    def fai   = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta.fai", checkIfExists: true)
    def target_bed = file("${launchDir}/tests/data/genomics/sarscov2/bed/test.bed", checkIfExists: true)
    input = [ [ id:'test'], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam.bai", checkIfExists: true) ]

    STRELKA_GERMLINE ( input, fasta, fai, target_bed )
}

