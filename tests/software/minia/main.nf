#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINIA } from '../../../software/minia/main.nf' addParams( options: [:] )

workflow test_minia {

    def input = []
    input = [ [ id:'test' ], // meta map
              [file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_2.fastq.gz", checkIfExists: true)] ]

    MINIA ( input )
}
