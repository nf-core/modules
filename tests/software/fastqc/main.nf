#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC as FASTQC_SE } from '../../../software/fastqc/main.nf'  addParams( options: [ publish_dir:'test_single_end' ] )
include { FASTQC as FASTQC_PE } from '../../../software/fastqc/main.nf'  addParams( options: [ publish_dir:'test_paired_end' ] )

/*
 * Test with single-end data
 */
workflow test_single_end {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/fastq/rna/test_single_end.fastq.gz", checkIfExists: true) ] ]

    FASTQC_SE ( input )
}

/*
 * Test with paired-end data
 */
workflow test_paired_end {

    def input = []
    input = [[id: 'test', single_end: false], // meta map
             [file("${launchDir}/tests/data/fastq/rna/test_R1.fastq.gz", checkIfExists: true),
              file("${launchDir}/tests/data/fastq/rna/test_R2.fastq.gz", checkIfExists: true)]]

    FASTQC_PE(input)

    emit:
    html = FASTQC_PE.out.html
    zip  = FASTQC_PE.out.zip

}
