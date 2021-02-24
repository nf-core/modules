#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQWISH_INDUCE } from '../../../../software/seqwish/induce/main.nf' addParams( options: [:] )

workflow test_seqwish_induce {

    def input = []
    input = [ [ id:'GCA_011545545.1_ASM1154554v1_cds_from_genomic' ], // meta map
              [ file("${launchDir}/tests/data/paf/GCA_011545545.1_ASM1154554v1_cds_from_genomic.paf", checkIfExists: true) ],
              [ file("${launchDir}/tests/data/fasta/sarscov2/GCA_011545545.1_ASM1154554v1_cds_from_genomic.fna", checkIfExists: true) ] ]

    SEQWISH_INDUCE ( input )
}