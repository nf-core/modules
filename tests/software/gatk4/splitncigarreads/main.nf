#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_SPLITNCIGARREADS } from '../../../../software/gatk4/splitncigarreads/main.nf' addParams( options: [:] )

workflow test_gatk4_splitncigarreads {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.bam", checkIfExists: true)] ]

    fasta = [ file("tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true),
              file("tests/data/genomics/sarscov2/fasta/test_genome.fasta.fai", checkIfExists: true),
              file("tests/data/genomics/sarscov2/fasta/test_genome.dict", checkIfExists: true)]

    GATK4_SPLITNCIGARREADS ( input, fasta )
}
