#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TIDDIT_SV }      from '../../../../software/tiddit/sv/main.nf'      addParams( options: [:] )

workflow test_tiddit_sv {
    def input = []
    def fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)
    def fai   = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta.fai", checkIfExists: true)

    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true) ] ]

    TIDDIT_SV ( input, fasta, fai )
}

workflow test_tiddit_sv_no_ref {
    def input = []

    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true) ] ]

    TIDDIT_SV ( input, [], [] )
}
