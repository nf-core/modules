#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { METHYLDACKEL_EXTRACT } from '../../../../software/methyldackel/extract/main.nf' addParams( options: [:] )

workflow test_methyldackel_extract {

    def input = []
    def fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genomic.fasta", checkIfExists: true)
    def fai   = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genomic.fasta.fai", checkIfExists: true)

    input     = [ [ id:'test', single_end:false ], // meta map
                  file("${launchDir}/tests/data/genomics/sarscov2/bam/test_methylated_paired_end.sorted.bam", checkIfExists: true),
                  file("${launchDir}/tests/data/genomics/sarscov2/bam/test_methylated_paired_end.sorted.bam.bai", checkIfExists: true) ]

    METHYLDACKEL_EXTRACT (
        input,
        fasta,
        fai
    )
}
