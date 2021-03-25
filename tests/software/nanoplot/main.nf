#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NANOPLOT  } from '../../../software/nanoplot/main.nf' addParams( options: [:] )

workflow test_nanoplot {
    def input = []
    input = [ [ id:'test' ], // meta map
            [ file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true) ],
            [ file(params.test_data['sarscov2']['nanopore']['test_sequencing_summary'], checkIfExists: true) ] ]
              
    NANOPLOT ( input )
}
