#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TIDDIT_SV }      from '../../../../software/tiddit/sv/main.nf'      addParams( options: [:] )
include { SAMTOOLS_FAIDX } from '../../../../software/samtools/faidx/main.nf' addParams( options: [:] )

workflow test_tiddit_sv {
    def input = []
    def fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)

    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true) ] ]

    SAMTOOLS_FAIDX ( fasta )

    TIDDIT_SV ( input, fasta, SAMTOOLS_FAIDX.out.fai )
}

workflow test_tiddit_sv_no_ref {
    def input = []
    def dummy_file  = file("${launchDir}/tests/data/dummy/dummy_file.txt", checkIfExists: true)
    def dummy_file2 = file("${launchDir}/tests/data/dummy/dummy_file2.txt", checkIfExists: true)

    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true) ] ]


    TIDDIT_SV ( input, dummy_file, dummy_file2 )
}