#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_COLLECTMULTIPLEMETRICS } from '../../../../software/picard/collectmultiplemetrics/main.nf' addParams( options: [:] )

workflow test_picard_collectmultiplemetrics {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/sarscov2_paired_aln.bam", checkIfExists: true) ]

    PICARD_COLLECTMULTIPLEMETRICS (
        input,
        file("${launchDir}/tests/data/genomics/sarscov2/fasta/GCA_011545545.1_ASM1154554v1_genomic.fasta", checkIfExists: true)
    )
}
