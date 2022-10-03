#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NANOPLOT  } from '../../../modules/nanoplot/main.nf'

workflow test_nanoplot_summary {
    def input = []
    input = [ [ id:'test' ], // meta map
            [ file(params.test_data['sarscov2']['nanopore']['test_sequencing_summary'], checkIfExists: true) ] ]
              
    NANOPLOT ( input )
}

workflow test_nanoplot_fastq {
    def input = []
    input = [ [ id:'test' ], // meta map
            [ file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true) ] ]

    NANOPLOT ( input )
}
