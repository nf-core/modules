#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { METHYLDACKEL_EXTRACT } from '../../../../software/methyldackel/extract/main.nf' addParams( options: [:] )

workflow test_methyldackel_extract {

    def input = []
    def fasta = file("${launchDir}/tests/data/generic/fasta/NC_010473.fa", checkIfExists: true)
    def fai   = file("${launchDir}/tests/data/generic/fasta/NC_010473.fa.fai", checkIfExists: true)

    input     = [ [ id:'test', single_end:false ], // meta map
                  file("${launchDir}/tests/data/generic/bam/test.paired_end_methylated.sorted.bam", checkIfExists: true),
                  file("${launchDir}/tests/data/generic/bam/test.paired_end_methylated.sorted.bam.bai", checkIfExists: true) ]

    METHYLDACKEL_EXTRACT ( input, fasta, fai )
}
