#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQWISH_INDUCE } from '../../../../software/seqwish/induce/main.nf' addParams( options: [:] )

workflow test_seqwish_induce {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/paf/test_cds_from_genomic.paf", checkIfExists: true) ],
              [ file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_cds_from_genomic.fasta", checkIfExists: true) ] ]

    SEQWISH_INDUCE ( input )
}
