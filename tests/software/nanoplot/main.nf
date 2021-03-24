#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NANOPLOT  } from '../../../software/nanoplot/main.nf' addParams( options: [:] )

workflow test_nanoplot {
    def input = []
    input = [ [ id:'test' ], // meta map
            [ file('tests/data/genomics/sarscov2/nanopore/fastq/test.fastq.gz', checkIfExists: true) ],
            [ file('tests/data/genomics/sarscov2/nanopore/sequencing_summary/test.sequencing_summary.txt', checkIfExists: true) ] ]
              
    NANOPLOT ( input )
}
