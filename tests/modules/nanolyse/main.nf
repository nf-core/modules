#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NANOLYSE } from '../../../modules/nanolyse/main.nf'

workflow test_nanolyse {
    input = [ 
        [ id:'test' ], // meta map
        [ file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true)]
    ]

    fasta = file("https://github.com/wdecoster/nanolyse/raw/master/reference/lambda.fasta.gz", checkIfExists: true)

    NANOLYSE ( input, fasta )
}
