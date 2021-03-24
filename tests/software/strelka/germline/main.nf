#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STRELKA_GERMLINE } from '../../../../software/strelka/germline/main.nf' addParams( options: [:] )

workflow test_strelka_germline {
    input   = [ [ id:'test'], // meta map
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam.bai", checkIfExists: true) 
              ]
    fasta   = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    fai     = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta.fai", checkIfExists: true)
    targets = []
    
    STRELKA_GERMLINE ( input, fasta, fai, targets )
}

workflow test_strelka_germline_target_bed {
    input   = [ [ id:'test'], // meta map
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam.bai", checkIfExists: true) 
              ]
    fasta   = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    fai     = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta.fai", checkIfExists: true)
    targets = file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed", checkIfExists: true)

    STRELKA_GERMLINE ( input, fasta, fai, targets )
}

