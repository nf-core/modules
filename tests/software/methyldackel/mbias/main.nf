#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { METHYLDACKEL_MBIAS } from '../../../../software/methyldackel/mbias/main.nf' addParams( options: [:] )

workflow test_methyldackel_mbias {

    def input = []
    def fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)
    def fai   = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta.fai", checkIfExists: true)

    input     = [ [ id:'test', single_end:false ], // meta map
                  file("${launchDir}/tests/data/genomics/sarscov2/bam/test_methylated_paired_end.sorted.bam", checkIfExists: true),
                  file("${launchDir}/tests/data/genomics/sarscov2/bam/test_methylated_paired_end.sorted.bam.bai", checkIfExists: true) ]

    METHYLDACKEL_MBIAS (
        input,
        fasta,
        fai
    )
}
